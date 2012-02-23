subroutine writewann
use modmain
use modldapu
use mod_nrkp
implicit none
integer i,j,ik,ikloc,n1,n2,j1,j2,nwan,ias
logical lpmat
complex(8), allocatable :: zm(:,:)
real(8), allocatable :: dm(:,:),eval(:)
call init0
call init1

wproc=mpi_grid_root()
if (.not.wannier) then
  write(*,*)
  write(*,'("Error(writewann_h) : WF generation is switched off")')
  write(*,*)
  call pstop
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
lpmat=.false.
if (task.eq.808) lpmat=.true.
call genwfnr(6,.false.,lmaxvr)
if (allocated(wann_h)) deallocate(wann_h)
allocate(wann_h(nwantot,nwantot,nkptnr))
wann_h=zzero
if (allocated(wann_e)) deallocate(wann_e)
allocate(wann_e(nwantot,nkptnr))
wann_e=0.d0
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  call genwann_h(.true.,evalsvnr(1,ik),wanncnrloc(1,1,ikloc),&
    &wann_h(1,1,ik),wann_e(1,ik))
enddo

!allocate(wann_ene_m(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
!allocate(wann_occ_m(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
!call wann_ene_occ_(wann_ene_m,wann_occ_m)

call mpi_grid_reduce(wann_h(1,1,1),nwantot*nwantot*nkptnr,dims=(/dim_k/),side=.true.)
call mpi_grid_reduce(wann_e(1,1),nwantot*nkptnr,dims=(/dim_k/),side=.true.)
!call mpi_grid_reduce(wann_p(1,1,1,1),3*nwantot*nwantot*nkpt,dims=(/dim_k/),side=.true.)
if (mpi_grid_root().and.task.eq.807) then
  call readfermi
  open(200,file="WANN_H.OUT",form="FORMATTED",status="REPLACE")
  write(200,'("# units of energy are Hartree, 1 Ha=",F18.10," eV")')ha2ev
  write(200,'("# fermi energy")')
  write(200,'(G18.10)')efermi
  write(200,'("# lattice vectors (3 rows)")')
  do i=1,3
    write(200,'(3G18.10)')avec(:,i)
  enddo
  write(200,'("# reciprocal lattice vectors (3 rows)")')
  do i=1,3
    write(200,'(3G18.10)')bvec(:,i)
  enddo
  write(200,'("# k-grid size")')
  write(200,'(3I6)')ngridk
  write(200,'("# number of k-points")')
  write(200,'(I8)')nkptnr
  write(200,'("# number of Wannier functions")')
  write(200,'(I8)')nwantot
  write(200,'("# number of atoms")')
  write(200,'(I8)')natmtot  
!  write(200,'("# occupancy matrix")')
!  do i=1,wann_natom
!    ias=wann_iprj(1,i)
!    write(200,'("#   atom : ",I4)')ias
!    do l=0,lmaxlu
!      if (sum(abs(wann_occ_m(idxlm(l,-l):idxlm(l,l),idxlm(l,-l):idxlm(l,l),:,:,ias))).gt.1d-8) then
!        t=0.0
!        do ispn=1,nspinor
!          write(200,'("#     ispn : ",I1)')ispn
!          write(200,'("#     real part")')
!          do lm1=l**2+1,(l+1)**2
!            write(200,'("#",1X,7F12.6)')(dreal(wann_occ_m(lm1,lm2,ispn,ispn,ias)),lm2=l**2+1,(l+1)**2)
!            t(ispn)=t(ispn)+dreal(wann_occ_m(lm1,lm1,ispn,ispn,ias))
!          enddo
!          write(200,'("#     imag part")')
!          do lm1=l**2+1,(l+1)**2
!            write(200,'("#",1X,7F12.6)')(dimag(wann_occ_m(lm1,lm2,ispn,ispn,ias)),lm2=l**2+1,(l+1)**2)
!          enddo
!          write(200,'("#     spin occupancy : ",F12.6)')t(ispn)
!        enddo !ispn
!        write(200,'("#   total occupancy : ",F12.6)')sum(t)
!      endif
!    enddo
!  enddo
  do ik=1,nkptnr
    write(200,'("# k-point : ",I8)')ik
    write(200,'("# weight")')
    write(200,'(G18.10)')wkptnr(ik)
    write(200,'("# lattice coordinates")')
    write(200,'(3G18.10)')vklnr(:,ik)
    write(200,'("# Cartesian coordinates")')
    write(200,'(3G18.10)')vkcnr(:,ik)
    write(200,'("# real part of H")')
    do i=1,nwantot
      write(200,'(255G18.10)')(dreal(wann_h(i,j,ik)),j=1,nwantot)
    enddo
    write(200,'("# imaginary part of H")')
    do i=1,nwantot
      write(200,'(255G18.10)')(dimag(wann_h(i,j,ik)),j=1,nwantot)
    enddo
    write(200,'("# eigen-values of H")')
    write(200,'(255G18.10)')(wann_e(j,ik),j=1,nwantot)
  enddo
  close(200)
endif

if (mpi_grid_root()) then
  do ias=1,natmtot
    nwan=nwannias(ias)
    if (nwan.ne.0) then
      write(*,*)"ias : ",ias,"  nwan : ",nwan
      allocate(zm(nwan,nwan))
      allocate(dm(nwan,nwan))
      allocate(eval(nwan))
      j1=0
      do n1=1,nwantot
        if (wan_info(wi_atom,n1).eq.ias) then
          j1=j1+1
          j2=0
          do n2=1,nwantot
            if (wan_info(wi_atom,n2).eq.ias) then
              j2=j2+1
              zm(j1,j2)=wann_h(n1,n2,1)
              dm(j1,j2)=dreal(zm(j1,j2))
            endif
          enddo
        endif
      enddo
      write(*,*)"Hamiltonian at gamma point :"
      do j1=1,nwan
        write(*,'(255F12.6)')(dreal(zm(j1,j2)),j2=1,nwan)
      enddo
      write(*,*)
      do j1=1,nwan
        write(*,'(255F12.6)')(dimag(zm(j1,j2)),j2=1,nwan)
      enddo
      write(*,*)"eigen-vectors : "
      call diagdsy(nwan,dm,eval)
      do j1=1,nwan
        write(*,'(2X,7G18.10)')(dm(j1,j2),j2=1,nwan)
      enddo
      write(*,*)
      write(*,'(2X,7G18.10)')(eval(j1),j1=1,nwan)
      write(*,*)
      write(*,*)"transpose of eigen-vectors : "
      do j1=1,nwan
        write(*,'(2X,7G18.10)')(dm(j2,j1),j2=1,nwan)
      enddo
      write(*,*)
      deallocate(zm,dm,eval)
    endif
  enddo
endif




!deallocate(wann_ene_m,wann_occ_m)
return
end
