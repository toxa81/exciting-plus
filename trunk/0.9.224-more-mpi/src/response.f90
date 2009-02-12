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

integer i,j,n,ngsh,gshmin,gshmax,ik,ikloc,ispn,istfv,ierr
character*100 fname
integer, external :: iknrglob
integer, external :: iknrglob2
logical, external :: in_set

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
  
if (task.eq.400) fname='RESPONSE_ME.OUT'
if (task.eq.401) fname='RESPONSE_CHI0.OUT'
if (task.eq.402) fname='RESPONSE_CHI.OUT'
if (task.eq.403.or.task.eq.404) fname='RESPONSE.OUT'

do_lr_io=.true.
#ifdef _MPI_
if (task.eq.403) do_lr_io=.false.
#endif

if (iproc.eq.0) then
  open(150,file=trim(fname),form='formatted',status='replace')
  write(150,'("Running on ",I4," proc.")')nproc
  if (.not.do_lr_io) write(150,'("I/O is turned off")')
endif

if (task.eq.400.or.task.eq.403) then
  allocate(igishell(ngvec))
  allocate(ishellng(ngvec,2))
  call getgshells(ngsh,igishell,ishellng)
! find minimum number of G-shells
  gshmin=100000
  gshmax=1
  do i=1,nvq0
! q-vector in lattice coordinates
    do j=1,3
      vq0l(j)=1.d0*ivq0m_list(j,i)/ngridk(j)+1d-12
    enddo
! find G-vector which brings q0 to first BZ
    vgq0l(:)=floor(vq0l(:))
    gshmin=min(gshmin,igishell(ivgig(vgq0l(1),vgq0l(2),vgq0l(3))))
    gshmax=max(gshmin,igishell(ivgig(vgq0l(1),vgq0l(2),vgq0l(3))))
  enddo !i
  if (iproc.eq.0) then
    write(150,*)
    write(150,'("Found minimum and maximum number of G-shells : ",2I4)')gshmin,gshmax
    write(150,'("Input minimum and maximum number of G-shells : ",2I4)')gshme1,gshme2
  endif
  if (gshmin.lt.gshme1) then
    if (iproc.eq.0) then
      write(150,*)
      write(150,'("Warning: minimum number of G-shells was changed from ",&
        &I4," to ",I4)')gshme1,gshmin
    endif
    gshme1=gshmin
  endif
  if (gshmax.gt.gshme2) then
    if (iproc.eq.0) then
      write(150,*)
      write(150,'("Warning: maximum number of G-shells was changed from ",&
        &I4," to ",I4)')gshme2,gshmax
    endif
    gshme2=gshmax
  endif
! test if G-shells are closed
  i=ishellng(gshme1,2)
  j=ishellng(gshme2,2)
  if (abs(gc(i)-gc(i+1)).lt.epslat.or.abs(gc(j)-gc(j+1)).lt.epslat) then
    write(*,*)
    write(*,'("Bug(response): G-shells are not closed")')
    write(*,*)
    call pstop
  endif
  if (gshme1.eq.1) then
    gvecme1=1
  else
    gvecme1=ishellng(gshme1-1,2)+1
  endif
  gvecme2=ishellng(gshme2,2)
  ngvecme=gvecme2-gvecme1+1
  if (iproc.eq.0) then
    write(150,*)
    write(150,'("G-shell limits      : ",2I4)')gshme1,gshme2
    write(150,'("G-vector limits     : ",2I4)')gvecme1,gvecme2
    write(150,'("number of G-vectors : ",I4)')ngvecme   
    call flushifc(150)
  endif
  deallocate(igishell)
  deallocate(ishellng)
endif

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

#ifdef _MPI_
if (iproc.eq.0) then
  write(150,*)
  write(150,'("Making MPI grid")')
endif
mpi_ndims=3
allocate(mpi_dims(mpi_ndims))
allocate(mpi_x(mpi_ndims))
allocate(mpi_periods(mpi_ndims))

if (nproc.le.nkptnr) then
  mpi_dims=(/nproc,1,1/)
else
  if (mod(nproc,nkptnr).ne.0.and.iproc.eq.0) then
    write(150,*)
    write(150,'("Warning: number of processors is not optimal")')
  endif
  mpi_dims=(/nkptnr,nproc/nkptnr,1/)
endif
mpi_dims=(/4,1,1/)
if (iproc.eq.0) then
  write(150,*)
  write(150,'("grid layout")')
  write(150,'("  k-points : ",I4)')mpi_dims(1)
  write(150,'("  G-vectors : ",I4)')mpi_dims(2)
  write(150,'("  q-points : ",I4)')mpi_dims(3)
  call flushifc(150)
endif
mpi_periods=.false.
call mpi_cart_create(MPI_COMM_WORLD,mpi_ndims,mpi_dims,mpi_periods,.false.,mpi_comm_cart,ierr)
call mpi_cart_get(mpi_comm_cart,mpi_ndims,mpi_dims,mpi_periods,mpi_x,ierr)
call mpi_cart_sub(mpi_comm_cart,(/.true.,.false.,.false./),mpi_comm_k,ierr)
call mpi_cart_sub(mpi_comm_cart,(/.false.,.true.,.true./),mpi_comm_gq,ierr)
call mpi_cart_sub(mpi_comm_cart,(/.false.,.true.,.false./),mpi_comm_g,ierr)
#endif

! distribute k-points over 1-st dimension of the grid
if (allocated(nkptnrloc)) deallocate(nkptnrloc)
allocate(nkptnrloc(0:mpi_dims(1)-1))
if (allocated(ikptnrloc)) deallocate(ikptnrloc)
allocate(ikptnrloc(0:mpi_dims(1)-1,2))
call splitk(nkptnr,mpi_dims(1),nkptnrloc,ikptnrloc)
nkptnr_loc=nkptnrloc(mpi_x(1))


!if (.true.) then
!  gvecme1=516
!  gvecme2=516
!  ngvecme=1
!  gvecchi1=516
!  gvecchi2=516
!  ngvecchi=1
!endif
  

if (task.eq.400.or.task.eq.401.or.task.eq.403.or.task.eq.404) then
! get occupancies and energies of states
  allocate(occsvnr(nstsv,nkptnr))
  allocate(evalsvnr(nstsv,nkptnr))
  call timer_start(3)
  if (iproc.eq.0) then
    write(150,*)
    write(150,'("Reading energies and occupancies of states")')
    call flushifc(150)
  endif
! if parallel I/O
#ifdef _PIO_
  occsvnr=0.d0
  evalsvnr=0.d0
! only subset of processors will read from file
  if (in_set((/1,0,0/))) then
    do ikloc=1,nkptnr_loc
      ik=iknrglob2(ikloc,mpi_x(1))
      call getoccsv(vklnr(1,ik),occsvnr(1,ik))
      call getevalsv(vklnr(1,ik),evalsvnr(1,ik))
    enddo
  endif
  call dsync2((/1,0,0/),occsvnr,nstsv*nkptnr,.true.,.true.)
  call dsync2((/1,0,0/),evalsvnr,nstsv*nkptnr,.true.,.true.)
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
  if (iproc.eq.0) then
    write(150,'("Done in ",F8.2," seconds")')timer(3,2)
    call flushifc(150)
  endif
endif

if (task.eq.400.or.task.eq.403) then
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
  if (iproc.eq.0) then
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
  if (in_set((/1,0,0/))) then
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
      call barrier2(mpi_comm_k)
    enddo
#endif
  endif !in_set
  call barrier
  call dsync2((/1,0,0/),wfsvmtloc, &
    lmmaxvr*nrfmax*natmtot*nstsv*nspinor*nkptnr_loc*2, &
    .false.,.false.)
  call dsync2((/1,0,0/),wfsvitloc, &
    ngkmax*nstsv*nspinor*nkptnr_loc*2,.false.,.false.)
  call timer_stop(1)
  if (iproc.eq.0) then
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
  do i=1,nvq0
    call response_me(ivq0m_list(1,i),wfsvmtloc,wfsvitloc,ngknr, &
      igkignr,occsvnr)
  enddo
endif

if (task.eq.401) then
! calculate chi0
  do i=1,nvq0
    call response_chi0(ivq0m_list(1,i),evalsvnr)
  enddo
endif

if (task.eq.402) then
! calculate chi
  do i=1,nvq0
    call response_chi(ivq0m_list(1,i))
  enddo
endif

if (task.eq.403) then
  do i=1,nvq0
    call response_me(ivq0m_list(1,i),wfsvmtloc,wfsvitloc,ngknr, &
      igkignr,occsvnr)
    call response_chi0(ivq0m_list(1,i),evalsvnr)
    call response_chi(ivq0m_list(1,i))
  enddo
endif

if (task.eq.404) call response_jdos(occsvnr,evalsvnr)

if (iproc.eq.0) close(150)

if (task.eq.400.or.task.eq.403) then
  deallocate(wfsvmtloc)
  deallocate(wfsvitloc)
  deallocate(ngknr)
  deallocate(igkignr)
  deallocate(occsvnr)
endif

return
end
