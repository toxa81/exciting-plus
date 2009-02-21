subroutine response
use modmain
#ifdef _MPI_
use mpi
#endif
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

complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: wfsvitloc(:,:,:,:)
complex(8), allocatable :: wfsvmtloc(:,:,:,:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
real(8), allocatable :: occsvnr(:,:)
real(8), allocatable :: evalsvnr(:,:)
complex(8), allocatable :: wfsvmt_t(:,:,:,:,:)
complex(8), allocatable :: wfc_t(:,:,:)

integer i,j,n,ngsh,gshmin,gshmax,ik,ikloc,ispn,istfv,ierr,rank
character*100 fname,qnm
integer, external :: iknrglob
integer, external :: iknrglob2
logical, external :: root_of

integer vgq0l(3)

! initialise universal variables
call init0
call init1

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

! make MPI cartesian grid
if (task.eq.400.or.task.eq.401) then
  mpi_ndims=3
endif
if (task.eq.402) then
  mpi_ndims=2
endif
if (allocated(mpi_dims)) deallocate(mpi_dims)
allocate(mpi_dims(mpi_ndims))
if (allocated(mpi_x)) deallocate(mpi_x)
allocate(mpi_x(mpi_ndims))
if (allocated(mpi_periods)) deallocate(mpi_periods)
allocate(mpi_periods(mpi_ndims))
mpi_periods=.false.
if (task.eq.400.or.task.eq.401) then
#ifdef _MPI_
  if (nproc.le.nkptnr) then
    mpi_dims=(/nproc,1,1/)
  else
    i=nproc/nkptnr
    if (i.le.nvq0) then
      mpi_dims=(/nkptnr,1,i/)
    else
      mpi_dims=(/nkptnr,nproc/(nkptnr*nvq0),nvq0/)
    endif
  endif
#else
  mpi_dims=(/1,1,1/)
#endif
endif

#ifdef _MPI_
call mpi_cart_create(MPI_COMM_WORLD,mpi_ndims,mpi_dims,mpi_periods,   &
  .false.,mpi_comm_cart,ierr)
call mpi_cart_get(mpi_comm_cart,mpi_ndims,mpi_dims,mpi_periods,mpi_x, &
  ierr)  
#endif

call qname(ivq0m_list(:,mpi_x(3)+1),qnm)
if (task.eq.400) then
  if (root_of((/0,0,1/))) then
    wproc=.true.
    fname=trim(qnm)//"_ME.OUT"
    open(150,file=trim(fname),form='formatted',status='replace')
  else
    wproc=.false.
  endif
endif
if (task.eq.401) then
  if (root_of((/0,0,1/))) then
    wproc=.true.
    fname=trim(qnm)//"_CHI0.OUT"
    open(150,file=trim(fname),form='formatted',status='replace')
  else
    wproc=.false.
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

#ifdef _MPI_
if (task.eq.400) then
  call mpi_cart_sub(mpi_comm_cart,(/.true.,.false.,.false./),mpi_comm_k,ierr)
  call mpi_cart_sub(mpi_comm_cart,(/.true.,.false.,.true./),mpi_comm_kq,ierr)
  call mpi_cart_sub(mpi_comm_cart,(/.false.,.true.,.true./),mpi_comm_gq,ierr)
  call mpi_cart_sub(mpi_comm_cart,(/.false.,.true.,.false./),mpi_comm_g,ierr)
endif
if (task.eq.401) then
  call mpi_cart_sub(mpi_comm_cart,(/.true.,.false.,.false./),mpi_comm_k,ierr)
  call mpi_cart_sub(mpi_comm_cart,(/.true.,.false.,.true./),mpi_comm_kq,ierr)
!  call mpi_cart_sub(mpi_comm_cart,(/.false.,.true.,.true./),mpi_comm_gq,ierr)
  call mpi_cart_sub(mpi_comm_cart,(/.false.,.true.,.false./),mpi_comm_b,ierr)
endif

#endif

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
  if (root_of((/1,0,0/))) then
    do ikloc=1,nkptnr_loc
      ik=iknrglob2(ikloc,mpi_x(1))
      call getoccsv(vklnr(1,ik),occsvnr(1,ik))
      call getevalsv(vklnr(1,ik),evalsvnr(1,ik))
    enddo
  endif
  call d_reduce_cart((/1,0,0/),.true.,evalsvnr,nstsv*nkptnr)
  call d_bcast((/0,1,1/),evalsvnr,nstsv*nkptnr)
  call d_reduce_cart((/1,0,0/),.true.,occsvnr,nstsv*nkptnr)
  call d_bcast((/0,1,1/),occsvnr,nstsv*nkptnr)
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
  allocate(evecfv(nmatmax,nstfv,nspnfv))
  allocate(evecsv(nstsv,nstsv))
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
  if (wproc) then
    write(150,*)
    write(150,'("Size of wfsvmt array (Mb) : ",I6)')                   &
      16*lmmaxvr*nrfmax*natmtot*nstsv*nspinor*nkptnrloc(0)/1024/1024
    write(150,'("Size of wfsvit array (Mb) : ",I6)')                   &
      16*ngkmax*nstsv*nspinor*nkptnrloc(0)/1024/1024
    write(150,*)
    write(150,'("Reading eigen-vectors")')
    call flushifc(150)
  endif
  call timer_start(1)
! read and transform eigen-vectors
  wfsvmtloc=0.d0
  wfsvitloc=0.d0
  if (root_of((/1,0,0/))) then
#ifndef _PIO_
    do i=0,mpi_dims(1)-1
      if (i.eq.mpi_x(1)) then
#endif
        do ikloc=1,nkptnr_loc
          ik=iknrglob2(ikloc,mpi_x(1))
          call getevecfv(vklnr(1,ik),vgklnr(1,1,ikloc),evecfv)
          call getevecsv(vklnr(1,ik),evecsv)
! get apw coeffs 
          call match(ngknr(ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc),         &
            sfacgknr(1,1,ikloc),apwalm)
! generate wave functions in muffin-tins
          call genwfsvmt(lmaxvr,lmmaxvr,ngknr(ikloc),evecfv,evecsv,apwalm, &
            wfsvmtloc(1,1,1,1,1,ikloc))
! generate wave functions in interstitial
          call genwfsvit(ngknr(ikloc),evecfv,evecsv,wfsvitloc(1,1,1,ikloc))
        enddo !ikloc
#ifndef _PIO_
      endif
      call grpbarrier(mpi_comm_k)
    enddo
#endif
  endif !root_of
  call barrier
  call d_bcast((/0,1,1/),wfsvmtloc, &
    lmmaxvr*nrfmax*natmtot*nstsv*nspinor*nkptnr_loc*2)
  call d_bcast((/0,1,1/),wfsvitloc, &
    ngkmax*nstsv*nspinor*nkptnr_loc*2)
  call timer_stop(1)
  if (wproc) then
    write(150,'("Done in ",F8.2," seconds")')timer(1,2)
    call flushifc(150)
  endif
  deallocate(evecfv,evecsv)
  deallocate(apwalm)
  deallocate(vgklnr)
  deallocate(vgkcnr)
  deallocate(gknr)
  deallocate(tpgknr)
  deallocate(sfacgknr)
endif

if (task.eq.400) then
! calculate matrix elements
  call response_me(ivq0m_list(1,mpi_x(3)+1),wfsvmtloc,wfsvitloc,ngknr, &
    igkignr,occsvnr)
endif

if (task.eq.401) then
! calculate chi0
  call response_chi0(ivq0m_list(1,mpi_x(3)+1),evalsvnr)
endif

if (task.eq.402) then
! calculate chi
  call response_chi(ivq0m_list(1,mpi_x(2)+1))
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
