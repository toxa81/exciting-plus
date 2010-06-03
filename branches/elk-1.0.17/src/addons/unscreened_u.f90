subroutine unscreened_u
use modmain
use mod_addons_q
use mod_nrkp
implicit none
real(8) t0,t1,t2
integer iq,ig,i,ik,ikloc,j
integer v1(3)
real(8) v2(3)
logical lgamma,wproc1
character*100 qnm

call init0
call init1
if (.not.mpi_grid_in()) return

if (.not.wannier) then
  write(*,*)
  write(*,'("Error(unscreened_u): Wannier functions are off")')
  write(*,*)
  call pstop
endif
wannier_megq=.true.
lgamma=.true.
call init_qbz(lgamma)
! create q-directories
if (mpi_grid_root()) then
  call system("mkdir -p q")
  do iq=1,nvq
    call getqdir(iq,vqm(:,iq),qnm)
    call system("mkdir -p "//trim(qnm))
  enddo
endif


if (allocated(vqlnr)) deallocate(vqlnr)
allocate(vqlnr(3,nvq))
if (allocated(vqcnr)) deallocate(vqcnr)
allocate(vqcnr(3,nvq))
if (allocated(vql)) deallocate(vql)
allocate(vql(3,nvq))
if (allocated(vqc)) deallocate(vqc)
allocate(vqc(3,nvq))
if (allocated(ig0q)) deallocate(ig0q)
allocate(ig0q(nvq))

t0=gqmax**2
ngqmax=0
do iq=1,nvq
  vqlnr(:,iq)=dble(vqm(:,iq))/ngridk(:)
  do ig=1,ngvec
    v1(:)=vqm(:,iq)-ngridk(:)*ivg(:,ig)
    if (v1(1).ge.0.and.v1(1).lt.ngridk(1).and.&
        v1(2).ge.0.and.v1(2).lt.ngridk(2).and.&
        v1(3).ge.0.and.v1(3).lt.ngridk(3)) then
      ig0q(iq)=ig
      vql(:,iq)=dble(v1(:))/ngridk(:)
      exit
    endif
  enddo !ig
  call r3mv(bvec,vqlnr(:,iq),vqcnr(:,iq))
  call r3mv(bvec,vql(:,iq),vqc(:,iq))
  i=0
  do ig=1,ngvec
    v2(:)=vgc(:,ig)+vqc(:,iq)
    t2=v2(1)**2+v2(2)**2+v2(3)**2
    if (t2.le.t0) i=i+1
  enddo
  ngqmax=max(ngqmax,i)
enddo
if (allocated(ngq)) deallocate(ngq)
allocate(ngq(nvq))
if (allocated(gq)) deallocate(gq)
allocate(gq(ngqmax,nvq))
if (allocated(vgqc)) deallocate(vgqc)
allocate(vgqc(3,ngqmax,nvq))
if (allocated(igqig)) deallocate(igqig)
allocate(igqig(ngqmax,nvq))
do iq=1,nvq
  ngq(iq)=0
  do ig=1,ngvec
    v2(:)=vgc(:,ig)+vqc(:,iq)
    t2=v2(1)**2+v2(2)**2+v2(3)**2
    if (t2.le.t0) then
      ngq(iq)=ngq(iq)+1
      gq(ngq(iq),iq)=sqrt(t2)
      vgqc(:,ngq(iq),iq)=v2(:)
      igqig(ngq(iq),iq)=ig
    endif
  enddo !ig
enddo    


wproc1=.false.
if (mpi_grid_root()) then
  wproc1=.true.
  open(151,file="UNSCREENED_U.OUT",form="FORMATTED",status="REPLACE")
  call timestamp(151)
endif
wproc=wproc1
if (wproc1) then
  write(151,'("Total number of q-vectors        : ",I6)')nvq
  write(151,'("Total number of processors       : ",I6)')nproc
  write(151,'("MPI grid size                    : ",3I6)')mpi_grid_size
  call flushifc(151)
endif

! read the density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
call getufr
call genufrp
! read Fermi energy
if (mpi_grid_root()) call readfermi
call mpi_grid_bcast(efermi)
! generate wave-functions for entire BZ
call genwfnr(151,.false.)
if (spinpol) then
  if (allocated(spinor_ud)) deallocate(spinor_ud)
  allocate(spinor_ud(2,nstsv,nkptnr))
  spinor_ud=0
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do j=1,nstsv
      t1=sum(abs(evecsvloc(1:nstfv,j,ikloc)))
      if (t1.gt.1d-10) spinor_ud(1,j,ik)=1
      t1=sum(abs(evecsvloc(nstfv+1:nstsv,j,ikloc)))
      if (t1.gt.1d-10) spinor_ud(2,j,ik)=1
    enddo
  enddo
  call mpi_grid_reduce(spinor_ud(1,1,1),2*nstsv*nkptnr,dims=(/dim_k/),all=.true.)
endif  
call getimegqwan(.true.)
if (wproc1) then
  write(151,*)
  write(151,'("Number of Wannier transitions : ",I6)')nmegqwan
  write(151,'("List of Wannier transitions (n n1 T) ")')
  do i=1,nmegqwan
    write(151,'(I4,4X,I4,4X,3I3)')imegqwan(:,i)
  enddo
  call timestamp(151)
  write(151,*)
  write(151,'("Translation limits : ",6I6)')megqwan_tlim(:,1), &
    megqwan_tlim(:,2),megqwan_tlim(:,3)
  call flushifc(151)
endif
if (wproc1) then
  close(151)
endif

return
end