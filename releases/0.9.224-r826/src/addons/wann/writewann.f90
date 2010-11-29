subroutine writewann
use modmain
use mod_mpi_grid
implicit none
integer i,j,ik,ikloc,idm,lm,ir,ispn,n,ias,io
complex(8), allocatable :: wann_rf(:,:,:,:)
real(8), allocatable :: eval(:)
character*20 fname

call init0
call init1

wproc=mpi_grid_root()
if (.not.wannier) then
  if (wproc) then
    write(*,*)
    write(*,'("Error(writewann_h) : WF generation is switched off")')
    write(*,*)
    call pstop
  endif
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

if (task.eq.601) then
  allocate(evecfvloc(nmatmax,nstfv,nspnfv,nkptloc))
  allocate(evecsvloc(nstsv,nstsv,nkptloc))
endif
if (task.eq.602) then
  allocate(wann_rf(lmmaxvr,nrmtmax,nspinor,nwann))
  wann_rf=zzero
endif
evalsv=0.d0
if (mpi_grid_side(dims=(/dim_k/))) then
  do i=0,mpi_grid_size(dim_k)-1
    if (mpi_grid_x(dim_k).eq.i) then
      do ikloc=1,nkptloc
        ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
        !call getwann(ikloc)
        if (task.eq.600) then
          call getevalsv(vkl(1,ik),evalsv(1,ik))
        endif
        if (task.eq.601) then
          call getevecfv(vkl(1,ik),vgkl(1,1,1,ikloc),evecfvloc(1,1,1,ikloc))
          call getevecsv(vkl(1,ik),evecsvloc(1,1,ikloc))
        endif
      enddo
    endif
    call mpi_grid_barrier(dims=(/dim_k/))
  enddo
endif

do ikloc=1,nkptloc
  if (task.eq.600) call genwann_h(ikloc)
  if (task.eq.601) call genwann_p(ikloc,evecfvloc(1,1,1,ikloc), &
    evecsvloc(1,1,ikloc))
  if (task.eq.602) then
    do n=1,nwann
      ias=iwann(1,n)
      do ispn=1,nspinor
        do ir=1,nrmt(ias2is(ias))
          do lm=1,lmmaxvr
            do io=1,nrfmax
              wann_rf(lm,ir,ispn,n)=wann_rf(lm,ir,ispn,n)+&
                urf(ir,lm2l(lm),io,ias)*wann_unkmt(lm,io,ias,ispn,n,ikloc)
            enddo
          enddo !lm
        enddo !ir
      enddo !ispn
    enddo !n
  endif
enddo
call mpi_grid_reduce(wann_h(1,1,1),nwann*nwann*nkpt,dims=(/dim_k/),side=.true.)
call mpi_grid_reduce(wann_p(1,1,1,1),3*nwann*nwann*nkpt,dims=(/dim_k/),side=.true.)
if (task.eq.602) then
  do n=1,nwann
    call mpi_grid_reduce(wann_rf(1,1,1,n),lmmaxvr*nrmtmax*nspinor,dims=(/dim_k/),side=.true.)
  enddo
endif
if (wproc.and.task.eq.600) then
  allocate(eval(nwann))
  call readfermi
  open(200,file='WANN_H.OUT',form='formatted',status='replace')
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
  write(200,'(I8)')nkpt
  write(200,'("# number of Wannier functions")')
  write(200,'(I8)')nwann
  do ik=1,nkpt
    write(200,'("# k-point : ",I8)')ik
    write(200,'("# weight")')
    write(200,'(G18.10)')wkpt(ik)
    write(200,'("# lattice coordinates")')
    write(200,'(3G18.10)')vkl(:,ik)
    write(200,'("# Cartesian coordinates")')
    write(200,'(3G18.10)')vkc(:,ik)
    write(200,'("# real part of H")')
    do i=1,nwann
      write(200,'(255G18.10)')(dreal(wann_h(i,j,ik)),j=1,nwann)
    enddo
    write(200,'("# imaginary part of H")')
    do i=1,nwann
      write(200,'(255G18.10)')(dimag(wann_h(i,j,ik)),j=1,nwann)
    enddo
    write(200,'("# eigen-values of H")')
    call diagzhe(nwann,wann_h(1,1,ik),eval)
    write(200,'(255G18.10)')(eval(j),j=1,nwann)
  enddo	
  close(200)
!  open(200,file='WANN_H0.OUT',form='formatted',status='replace')
!  do i=1,nwann
!    write(200,'(6X,255G18.10)')(dreal(sum(wann_h(i,j,:))/nkpt),j=1,nwann)
!  enddo
!  write(200,*)
!  do i=1,nwann
!    write(200,'(6X,255G18.10)')(dimag(sum(wann_h(i,j,:))/nkpt),j=1,nwann)
!  enddo
!  close(200)
!  do i=1,nwann
!    wann_h(i,i,:)=wann_h(i,i,:)-efermi
!  enddo
!  wann_h=wann_h*ha2ev
!  open(200,file='WANN_H_OLD.OUT',form='formatted',status='replace')
!  write(200,*)nkpt,nwann
!  do ik=1,nkpt
!    write(200,*)1.d0 !wtkp(ikp)
!    do i=1,nwann
!      do j=1,nwann
!        write(200,*)dreal(wann_h(i,j,ik)),dimag(wann_h(i,j,ik))
!      enddo
!    enddo
!  enddo	
!  close(200)
  deallocate(eval)
endif
if (wproc.and.task.eq.601) then
  open(200,file='WANN_P.OUT',form='formatted',status='replace')
  do ik=1,nkpt
    write(200,'("ik : ",I4)')ik
    do idm=1,3
      write(200,'("  x : ",I4)')idm
      do i=1,nwann
        write(200,'(6X,255G18.10)')(dreal(wann_p(idm,i,j,ik)),j=1,nwann)
      enddo
      write(200,*)
      do i=1,nwann
        write(200,'(6X,255G18.10)')(dimag(wann_p(idm,i,j,ik)),j=1,nwann)
      enddo
    enddo
  enddo
  write(200,*)
  do idm=1,3
    write(200,'("  x : ",I4)')idm
    do i=1,nwann
      write(200,'(6X,255G18.10)')(dreal(sum(wann_p(idm,i,j,:))/nkpt),j=1,nwann)
    enddo
    write(200,*)
    do i=1,nwann
      write(200,'(6X,255G18.10)')(dimag(sum(wann_p(idm,i,j,:))/nkpt),j=1,nwann)
    enddo
  enddo
  close(200)
endif
if (wproc.and.task.eq.602) then
  do n=1,nwann
    write(fname,'("WANN_",I3.3,"_rfmt.OUT")')n
    open(200,file=trim(fname),form='formatted',status='replace')
    do lm=1,16
      do ir=1,nrmt(ias2is(iwann(1,n)))
        write(200,'(2G18.10)')spr(ir,ias2is(iwann(1,n))),abs(wann_rf(lm,ir,1,n))
      enddo
      write(200,*)
    enddo
    close(200)
  enddo
endif
if (task.eq.601) then
  deallocate(evecfvloc)
  deallocate(evecsvloc)
endif
if (task.eq.602) then
  deallocate(wann_rf)
endif
return
end
