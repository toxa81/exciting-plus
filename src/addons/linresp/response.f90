#ifdef _HDF5_
subroutine response
use modmain
use hdf5
use mod_nrkp
implicit none
#ifdef _PAPI_
integer ierr
include 'f90papi.h'
real real_time,cpu_time,mflops
integer*8 fp_ins
#endif

integer i,j,n,ik,ikloc,iq
integer nwloc
integer nvq0loc,iqloc,i1,i2,i3
character*100 qnm
real(8) t1
logical lgamma,wproc1,lpmat
logical, external :: wann_diel

! Task list:
!   400 - response functions
!   401 - cRPA
!   402 - bare U
!
! MPI grid:
!   (matrix elements) : (1) k-points x (2) interband transitions x 
!                     x (3) q-points 
!   (response) : (1) energy mesh x (2) number of fxc kernels x (3) q-points
!
! TODO: for U(w) change number of G-vectors from G-shells to |G+q| cutoff  
!
! TODO: variable q0 is confusing -> change to q
!
if (lrtype.eq.1.and..not.spinpol) then
  write(*,*)
  write(*,'("Error(response): can''t do magnetic response for unpolarized ground state")')
  write(*,*)
  call pstop
endif
if (nvq0.eq.0.and..not.(task.eq.401.or.task.eq.402)) then
  write(*,*)
  write(*,'("Error(response): no q-vectors")')
  write(*,*)
  call pstop
endif   
if (.not.wannier) wannier_chi0_chi=.false.
if (.not.spinpol) megqwan_afm=.false.
if (wannier_chi0_chi.or.task.eq.401.or.task.eq.402) wannier_megq=.true.

! this is enough for matrix elements
!lmaxvr=5

if (iproc.eq.0) call timestamp(6,'before init0')
! initialise universal variables
call init0
call init1
call mpi_world_barrier
if (iproc.eq.0) call timestamp(6,'after init1')
if (.not.mpi_grid_in()) return

! check if all q-vectors in BZ are required 
lgamma=.true.
if (task.eq.401.or.task.eq.402) then
  call init_qbz(lgamma)
  nfxca=1
  fxca0=0.d0
  fxca1=0.d0
endif

! check if momentum matrix is required
lpmat=.false.
do j=1,nvq0
  if (ivq0m_list(1,j).eq.0.and.ivq0m_list(2,j).eq.0.and. &
      ivq0m_list(3,j).eq.0) then
    lpmat=.true.
  endif
enddo
! create q-directories
if (mpi_grid_root()) then
  call system("mkdir -p q")
  do iq=1,nvq0
    call getqdir(iq,ivq0m_list(:,iq),qnm)
    call system("mkdir -p "//trim(qnm))
  enddo
endif

wproc1=.false.
if (mpi_grid_root()) then
  wproc1=.true.
  open(151,file="RESPONSE.OUT",form="FORMATTED",status="REPLACE")
  call timestamp(151)
endif
wproc=wproc1
if (wproc1) then
  write(151,'("Total number of q-vectors        : ",I6)')nvq0
  write(151,'("Total number of processors       : ",I6)')nproc
  write(151,'("MPI grid size                    : ",3I6)')mpi_grid_size
  if (nproc.gt.1) then
    write(151,'("Parallel file reading            : ",L1)')parallel_read
  endif
  write(151,'("Wannier functions                : ",L1)')wannier
  write(151,'("Matrix elements in Wannier basis : ",L1)')wannier_megq
  write(151,'("Response in Wannier basis        : ",L1)')wannier_chi0_chi
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
call geturf
call genurfprod
! read Fermi energy
if (mpi_grid_root()) call readfermi
call mpi_grid_bcast(efermi)
! generate wave-functions for entire BZ
call genwfnr(151,lpmat)
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

if (wannier_megq) then
  all_wan_ibt=(task.eq.401.or.task.eq.402)
  call getimegqwan(all_wan_ibt)
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
endif

! setup energy mesh
lr_dw=(lr_w1-lr_w0)/(lr_nw-1)
if (allocated(lr_w)) deallocate(lr_w)
allocate(lr_w(lr_nw))
do i=1,lr_nw
  lr_w(i)=dcmplx(lr_w0+lr_dw*(i-1),lr_eta)/ha2ev
enddo

if (task.eq.401.or.task.eq.402) then
  maxtr_uscrn=1
  ntr_uscrn=(2*maxtr_uscrn+1)**3
  if (allocated(vtl_uscrn)) deallocate(vtl_uscrn)
  allocate(vtl_uscrn(3,ntr_uscrn))
  if (allocated(ivtit_uscrn)) deallocate(ivtit_uscrn)
  allocate(ivtit_uscrn(-maxtr_uscrn:maxtr_uscrn,-maxtr_uscrn:maxtr_uscrn,&
    -maxtr_uscrn:maxtr_uscrn))
  n=0
  do i1=-maxtr_uscrn,maxtr_uscrn
    do i2=-maxtr_uscrn,maxtr_uscrn
      do i3=-maxtr_uscrn,maxtr_uscrn
        n=n+1
        vtl_uscrn(:,n)=(/i1,i2,i3/)
        ivtit_uscrn(i1,i2,i3)=n
      enddo
    enddo
  enddo
  nwloc=mpi_grid_map(lr_nw,dim_k)
  if (allocated(uscrnwan)) deallocate(uscrnwan)
  allocate(uscrnwan(nwann,nwann,ntr_uscrn,nwloc))
  uscrnwan=zzero
  if (allocated(ubarewan)) deallocate(ubarewan)
  allocate(ubarewan(nwann,nwann,ntr_uscrn))
  ubarewan=zzero
endif

! distribute q-vectors along 3-rd dimention
nvq0loc=mpi_grid_map(nvq0,dim_q)

#ifdef _PAPI_
call PAPIF_flops(real_time,cpu_time,fp_ins,mflops,ierr)
#endif

! main loop over q-points
do iqloc=1,nvq0loc
  iq=mpi_grid_map(nvq0,dim_q,loc=iqloc)
  call genmegq(iq,.true.)
  if (task.eq.400.or.task.eq.401) call genchi(iq)
  if (task.eq.402) call genubare(iq)
enddo
if (task.eq.401) call write_uscrn
if (task.eq.402) call write_ubare

#ifdef _PAPI_
call PAPIF_flops(real_time,cpu_time,fp_ins,mflops,ierr)
t1=dble(mflops)
call mpi_grid_reduce(t1)
#endif

if (wproc1) then
  write(151,*)
#ifdef _PAPI_
  write(151,'("Average performance (Gflops/proc) : ",F15.3)')t1/mpi_grid_nproc/1000.d0
#endif  
  write(151,'("Done.")')
  close(151)
endif

deallocate(wfsvmtloc)
deallocate(wfsvitloc)
deallocate(evecfvloc)
deallocate(evecsvloc)
deallocate(ngknr)
deallocate(igkignr)
deallocate(occsvnr)
deallocate(evalsvnr)   
if (wannier_megq) deallocate(wann_c)
return
end
#endif
