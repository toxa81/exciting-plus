subroutine response
use modmain
#ifdef _MPI_
use mpi
#endif
use hdf5
implicit none

integer, allocatable :: igishell(:)
integer, allocatable :: ishellng(:,:)

integer, allocatable :: ngknr(:)
integer, allocatable :: igkignr(:,:)
real(8), allocatable :: vgklnr(:,:,:)
real(8), allocatable :: vgkcnr(:,:)
real(8), allocatable :: gknr(:,:)
real(8), allocatable :: tpgknr(:,:,:)
complex(8), allocatable :: sfacgknr(:,:,:)

!complex(8), allocatable :: evecfv(:,:,:)
!complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: wfsvitloc(:,:,:,:)
complex(8), allocatable :: wfsvmtloc(:,:,:,:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
real(8), allocatable :: occsvnr(:,:)
real(8), allocatable :: evalsvnr(:,:)
complex(8), allocatable :: wfsvmt_t(:,:,:,:,:)
complex(8), allocatable :: wfc_t(:,:,:)

integer i,j,n,ngsh,gshmin,gshmax,ik,ikloc,ispn,istfv,ierr,rank
integer sz
character*100 fname,qnm
character*2 c2
integer, external :: iknrglob
integer, external :: iknrglob2
logical, external :: root_cart

integer vgq0l(3)

! initialise universal variables
call init0
call init1
! initialize HDF5 library
call h5open_f(ierr)

if (ncmag) then
  write(*,*)
  write(*,'("Error(response): can''t do response + non-collinear magnetism")')
  write(*,*)
  call pstop
endif

if (.not.spinpol) then
  spin_me=1
endif

if (lrtype.eq.1.and..not.spinpol) then
  write(*,*)
  write(*,'("Error(response): can''t do magnetic response for unpolarized ground state")')
  write(*,*)
  call pstop
endif

do_lr_io=.true.

call qname(ivq0m_list(:,mpi_x(3)+1),qnm)
wproc=.false.
if (task.eq.400) then
  if (root_cart((/1,1,0/))) then
    wproc=.true.
    fname=trim(qnm)//"_ME.OUT"
    open(150,file=trim(fname),form='formatted',status='replace')
  endif
endif
if (task.eq.401) then
  if (root_cart((/1,1,0/))) then
    wproc=.true.
    fname=trim(qnm)//"_CHI0.OUT"
    open(150,file=trim(fname),form='formatted',status='replace')
  endif
endif
if (task.eq.402) then
  write(c2,'(I2.2)')mpi_x(2)
  if (root_cart((/1,0,0/))) then
    wproc=.true.
    fname=trim(qnm)//"_CHI_"//c2//".OUT"
    open(150,file=trim(fname),form='formatted',status='replace')
  endif
endif

if (wproc) then
  write(150,'("Running on ",I4," proc.")')nproc
  if (.not.do_lr_io) write(150,'("Output of intermediate files is turned off")')
#ifdef _PIO_
  if (nproc.gt.1) then
    write(150,'("Reading files in parallel")')
  endif
#endif
  write(150,'("MPI grid size : ",3I4)')mpi_dims
endif

! distribute k-points over 1-st dimension of the grid
if (allocated(nkptnrloc)) deallocate(nkptnrloc)
allocate(nkptnrloc(0:mpi_dims(1)-1))
if (allocated(ikptnrloc)) deallocate(ikptnrloc)
allocate(ikptnrloc(0:mpi_dims(1)-1,2))
call splitk(nkptnr,mpi_dims(1),nkptnrloc,ikptnrloc)
nkptnr_loc=nkptnrloc(mpi_x(1))

if (task.eq.402.or.task.eq.403) then
  allocate(igishell(ngvec))
  allocate(ishellng(ngvec,2))
  call getgshells(ngsh,igishell,ishellng)
  if (gshchi1.eq.1) then
    gvecchi1=1
  else
    gvecchi1=ishellng(gshchi1-1,2)+1
  endif
  gvecchi2=ishellng(gshchi2,2)
  deallocate(igishell)
  deallocate(ishellng)
endif



!if (.true.) then
!  gvecme1=516
!  gvecme2=516
!  ngvecme=1
!  gvecchi1=516
!  gvecchi2=516
!  ngvecchi=1
!endif
  

if (task.eq.400.or.task.eq.401.or.task.eq.404) then
! get occupancies and energies of states
  allocate(occsvnr(nstsv,nkptnr))
  allocate(evalsvnr(nstsv,nkptnr))
  call timer_start(3)
  if (wproc) then
    write(150,*)
    write(150,'("Reading energies and occupancies of states")')
    call flushifc(150)
  endif
! if parallel I/O
#ifdef _PIO_
  occsvnr=0.d0
  evalsvnr=0.d0
! only subset of processors will read from file
  if (root_cart((/0,1,1/))) then
    do ikloc=1,nkptnr_loc
      ik=iknrglob2(ikloc,mpi_x(1))
      call getoccsv(vklnr(1,ik),occsvnr(1,ik))
      call getevalsv(vklnr(1,ik),evalsvnr(1,ik))
    enddo
    call d_reduce_cart(comm_cart_100,.true.,evalsvnr,nstsv*nkptnr)
    call d_reduce_cart(comm_cart_100,.true.,occsvnr,nstsv*nkptnr)
  endif
  call d_bcast_cart(comm_cart_011,evalsvnr,nstsv*nkptnr)
  call d_bcast_cart(comm_cart_011,occsvnr,nstsv*nkptnr)
! if not parallel I/O
#else
  if (iproc.eq.0) then 
    do ik=1,nkptnr
      call getoccsv(vklnr(1,ik),occsvnr(1,ik))
      call getevalsv(vklnr(1,ik),evalsvnr(1,ik))
    enddo
  endif
  call dsync(occsvnr,nstsv*nkptnr,.false.,.true.)
  call dsync(evalsvnr,nstsv*nkptnr,.false.,.true.)
#endif
  call timer_stop(3)
  if (wproc) then
    write(150,'("Done in ",F8.2," seconds")')timer(3,2)
    call flushifc(150)
  endif
endif

if (task.eq.400) then
! read the density and potentials from file
  call readstate
! read Fermi energy from file
  call readfermi
! find the new linearisation energies
  call linengy
! generate the APW radial functions
  call genapwfr
! generate the local-orbital radial functions
  call genlofr
  call geturf
  call genurfprod
! generate G+k vectors for entire BZ (this is required to compute 
!   wave-functions at each k-point)
  allocate(vgklnr(3,ngkmax,nkptnr_loc))
  allocate(vgkcnr(3,ngkmax))
  allocate(gknr(ngkmax,nkptnr_loc))
  allocate(tpgknr(2,ngkmax,nkptnr_loc))
  allocate(ngknr(nkptnr_loc))
  allocate(sfacgknr(ngkmax,natmtot,nkptnr_loc))
  allocate(igkignr(ngkmax,nkptnr_loc))
  do ikloc=1,nkptnr_loc
    ik=iknrglob2(ikloc,mpi_x(1))
    call gengpvec(vklnr(1,ik),vkcnr(1,ik),ngknr(ikloc),igkignr(1,ikloc), &
      vgklnr(1,1,ikloc),vgkcnr,gknr(1,ikloc),tpgknr(1,1,ikloc))
    call gensfacgp(ngknr(ikloc),vgkcnr,ngkmax,sfacgknr(1,1,ikloc))
  enddo
  allocate(wfsvmtloc(lmmaxvr,nrfmax,natmtot,nstsv,nspinor,nkptnr_loc))
  allocate(wfsvitloc(ngkmax,nstsv,nspinor,nkptnr_loc))
  allocate(evecfvloc(nmatmax,nstfv,nspnfv,nkptnr_loc))
  allocate(evecsvloc(nstsv,nstsv,nkptnr_loc))
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
  if (wproc) then
    sz=lmmaxvr*nrfmax*natmtot*nstsv*nspinor
    sz=sz+ngkmax*nstsv*nspinor
    sz=sz+nmatmax*nstfv*nspnfv
    sz=sz+nstsv*nstsv
    sz=16*sz*nkptnr_loc/1024/1024
    write(150,*)
    write(150,'("Size of wave-function arrays (MB) : ",I6)')sz
    write(150,*)
    write(150,'("Reading eigen-vectors")')
    call flushifc(150)
  endif
  call timer_start(1)
! read and transform eigen-vectors
  wfsvmtloc=0.d0
  wfsvitloc=0.d0
  if (root_cart((/0,1,1/))) then
#ifndef _PIO_
    do i=0,mpi_dims(1)-1
      if (i.eq.mpi_x(1)) then
#endif
        do ikloc=1,nkptnr_loc
          ik=iknrglob2(ikloc,mpi_x(1))
          call getevecfv(vklnr(1,ik),vgklnr(1,1,ikloc),evecfvloc(1,1,1,ikloc))
          call getevecsv(vklnr(1,ik),evecsvloc(1,1,ikloc))
! get apw coeffs 
          call match(ngknr(ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc),        &
            sfacgknr(1,1,ikloc),apwalm)
! generate wave functions in muffin-tins
          call genwfsvmt(lmaxvr,lmmaxvr,ngknr(ikloc),evecfvloc(1,1,1,ikloc), &
            evecsvloc(1,1,ikloc),apwalm,wfsvmtloc(1,1,1,1,1,ikloc))
! generate wave functions in interstitial
          call genwfsvit(ngknr(ikloc),evecfvloc(1,1,1,ikloc), &
            evecsvloc(1,1,ikloc),wfsvitloc(1,1,1,ikloc))
        enddo !ikloc
#ifndef _PIO_
      endif
      call barrier(comm_cart_100)
    enddo
#endif
  endif !root_cart
  call barrier(comm_world)
  call d_bcast_cart(comm_cart_011,wfsvmtloc, &
    lmmaxvr*nrfmax*natmtot*nstsv*nspinor*nkptnr_loc*2)
  call d_bcast_cart(comm_cart_011,wfsvitloc, &
    ngkmax*nstsv*nspinor*nkptnr_loc*2)
  call timer_stop(1)
  if (wproc) then
    write(150,'("Done in ",F8.2," seconds")')timer(1,2)
    call flushifc(150)
  endif
!  deallocate(evecfv,evecsv)
  deallocate(apwalm)
  deallocate(vgklnr)
  deallocate(vgkcnr)
  deallocate(gknr)
  deallocate(tpgknr)
  deallocate(sfacgknr)
endif

if (task.eq.400) then
! calculate matrix elements
  call timer_start(10)
  call response_me(ivq0m_list(1,mpi_x(3)+1),wfsvmtloc,wfsvitloc,ngknr, &
    igkignr,occsvnr)
  call timer_stop(10)
  if (wproc) then
    write(150,'("Total time : ",F8.2," seconds")')timer(10,2)
    call flushifc(150)
  endif
endif

if (task.eq.401) then
! calculate chi0
  call timer_start(11)
  call response_chi0(ivq0m_list(1,mpi_x(3)+1),evalsvnr)
  call timer_stop(11)
  if (wproc) then
    write(150,'("Total time : ",F8.2," seconds")')timer(11,2)
    call flushifc(150)
  endif
endif

if (task.eq.402) then
! calculate chi
  call timer_start(12)
  call response_chi(ivq0m_list(1,mpi_x(3)+1))
  call timer_stop(12)
  if (wproc) then
    write(150,'("Total time : ",F8.2," seconds")')timer(12,2)
    call flushifc(150)
  endif
endif

if (task.eq.404) call response_jdos(occsvnr,evalsvnr)

if (wproc) close(150)

if (task.eq.400.or.task.eq.403) then
  deallocate(wfsvmtloc)
  deallocate(wfsvitloc)
  deallocate(ngknr)
  deallocate(igkignr)
  deallocate(occsvnr)
endif

return
end
