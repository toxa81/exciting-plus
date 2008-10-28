subroutine response
use modmain
implicit none

integer, allocatable :: igishell(:)
integer, allocatable :: ishellng(:,:)

integer, allocatable :: nkptlocnr(:)
integer, allocatable :: ikptlocnr(:,:)
integer, allocatable :: ngknr(:)
integer, allocatable :: igkignr(:,:)
real(8), allocatable :: vgklnr(:,:,:)
real(8), allocatable :: vgkcnr(:,:)
real(8), allocatable :: gknr(:,:)
real(8), allocatable :: tpgknr(:,:,:)
complex(8), allocatable :: sfacgknr(:,:,:)

complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: wfsvitloc(:,:,:)
complex(8), allocatable :: wfsvmtloc(:,:,:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
real(8), allocatable :: occsvnr(:,:)
real(8), allocatable :: wfnrmdev(:,:)

integer i,j,ngsh,gshmin,gshmax,gvecme1,gvecme2,ngvecme,gvecchi1,gvecchi2,ngvecchi,ik,ikloc,ig, &
  ispn,istfv
complex(8) zt1
character*100 fname

real(8) vq0l(3)
integer vgq0l(3)

! initialise universal variables
call init0
call init1

if (ndmag.eq.3) then
  write(*,*)
  write(*,'("Error(response): can''t do response + non-collinear magnetism")')
  write(*,*)
  call pstop
endif

if (.not.spinpol) then
  spin_me=1
  spin_chi=1
endif

if (task.eq.400) fname='RESPONSE_ME.OUT'
if (task.eq.401) fname='RESPONSE_CHI0.OUT'
if (task.eq.402) fname='RESPONSE_CHI.OUT'
if (task.eq.403) fname='RESPONSE.OUT'
if (task.eq.404) fname='RESPONSE.OUT'

if (iproc.eq.0) then
  open(150,file=trim(fname),form='formatted',status='replace')
  if (ismpi) then
    write(150,'("Running in parallel mode on ",I4," proc.")')nproc
  else
    write(150,'("Running in serial mode")')
  endif
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

allocate(nkptlocnr(0:nproc-1))
allocate(ikptlocnr(0:nproc-1,2))
call splitk(nkptnr,nproc,nkptlocnr,ikptlocnr)

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
  allocate(vgklnr(3,ngkmax,nkptlocnr(iproc)))
  allocate(vgkcnr(3,ngkmax))
  allocate(gknr(ngkmax,nkptlocnr(iproc)))
  allocate(tpgknr(2,ngkmax,nkptlocnr(iproc)))
  allocate(ngknr(nkptlocnr(iproc)))
  allocate(sfacgknr(ngkmax,natmtot,nkptlocnr(iproc)))
  allocate(igkignr(ngkmax,nkptlocnr(iproc)))
  do ikloc=1,nkptlocnr(iproc)
    ik=ikptlocnr(iproc,1)+ikloc-1
    call gengpvec(vklnr(1,ik),vkcnr(1,ik),ngknr(ikloc),igkignr(1,ikloc), &
      vgklnr(1,1,ikloc),vgkcnr,gknr(1,ikloc),tpgknr(1,1,ikloc))
    call gensfacgp(ngknr(ikloc),vgkcnr,ngkmax,sfacgknr(1,1,ikloc))
  enddo
  allocate(wfsvmtloc(lmmaxvr,nrfmax,natmtot,nstsv,nkptlocnr(iproc)))
  allocate(wfsvitloc(nmatmax,nstsv,nkptlocnr(iproc)))
  allocate(evecfv(nmatmax,nstfv,nspnfv))
  allocate(evecsv(nstsv,nstsv))
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
  if (iproc.eq.0) then
    write(150,*)
    write(150,'("Size of wfsvmt array (Mb) : ",I6)')lmmaxvr*nrfmax*natmtot*nstsv*nkptlocnr(0)/1024/1024
    write(150,'("Size of wfsvit array (Mb) : ",I6)')nmatmax*nstsv*nkptlocnr(0)/1024/1024
    write(150,*)
    write(150,'("Reading eigen-vectors")')
    call flushifc(150)
  endif
  allocate(wfnrmdev(nstsv*(nstsv+1)/2,nkptnr))
  wfnrmdev=0.d0
! read and transform eigen-vectors
  do ikloc=1,nkptlocnr(0)
    do i=0,nproc-1
      if (iproc.eq.i.and.ikloc.le.nkptlocnr(iproc)) then
        ik=ikptlocnr(iproc,1)+ikloc-1
        call getevecfv(vklnr(1,ik),vgklnr(1,1,ikloc),evecfv)
        call getevecsv(vklnr(1,ik),evecsv)
      endif
      call barrier
    enddo !i
    if (ikloc.le.nkptlocnr(iproc)) then
! get apw coeffs 
      call match(ngknr(ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc), &
        sfacgknr(1,1,ikloc),apwalm)
      call genwfsvmt(lmaxvr,lmmaxvr,ngknr(ikloc),evecfv,evecsv,apwalm, &
        wfsvmtloc(1,1,1,1,ikloc))
      call genwfsvit(ngknr(ikloc),evecfv,evecsv,wfsvitloc(1,1,ikloc))
      call wfsvprodk(ngknr(ikloc),igkignr(1,ikloc),wfsvmtloc(1,1,1,1,ikloc), &
        wfsvitloc(1,1,ikloc),wfnrmdev(1,ik))
    endif
  enddo !ikloc
  call dsync(wfnrmdev,nkptnr*(nstsv*(nstsv+1)/2),.true.,.false.)
  if (iproc.eq.0) then
    write(150,'("Done.")')
    write(150,*)
    write(150,'("Maximum WF norm deviation : ",G18.10)')maxval(wfnrmdev)
    write(150,'("Average WF norm deviation : ",G18.10)')sum(wfnrmdev)/nkptnr/(nstsv*(nstsv+1)/2.d0)
    call flushifc(150)
  endif
  deallocate(wfnrmdev)
  deallocate(evecfv,evecsv)
  deallocate(apwalm)
  deallocate(vgklnr)
  deallocate(vgkcnr)
  deallocate(gknr)
  deallocate(tpgknr)
  deallocate(sfacgknr)
! get occupancy of states
  allocate(occsvnr(nstsv,nkptnr))
  if (iproc.eq.0) then 
    do ik=1,nkptnr
      call getoccsv(vklnr(1,ik),occsvnr(1,ik))
    enddo
  endif
  call dsync(occsvnr,nstsv*nkptnr,.false.,.true.)
endif

if (task.eq.400) then
! calculate matrix elements
  do i=1,nvq0
    call response_me(ivq0m_list(1,i),gvecme1,gvecme2,ngvecme,ikptlocnr, &
      nkptlocnr,wfsvmtloc,wfsvitloc,ngknr,igkignr,occsvnr)
  enddo
endif

if (task.eq.401) then
! calculate chi0
  do i=1,nvq0
    call response_chi0(ivq0m_list(1,i),ikptlocnr,nkptlocnr)
  enddo
endif

if (task.eq.402) then
! calculate chi
  do i=1,nvq0
    call response_chi(ivq0m_list(1,i),gvecchi1,gvecchi2)
  enddo
endif

if (task.eq.403) then
  do i=1,nvq0
    call response_me(ivq0m_list(1,i),gvecme1,gvecme2,ngvecme,ikptlocnr, &
      nkptlocnr,wfsvmtloc,wfsvitloc,ngknr,igkignr,occsvnr)
    call response_chi0(ivq0m_list(1,i),ikptlocnr,nkptlocnr)
    call response_chi(ivq0m_list(1,i),gvecchi1,gvecchi2)
  enddo
endif

if (task.eq.404) then
  call response_jdos
endif

if (iproc.eq.0) close(150)

deallocate(nkptlocnr)
deallocate(ikptlocnr)
if (task.eq.400.or.task.eq.403) then
  deallocate(wfsvmtloc)
  deallocate(wfsvitloc)
  deallocate(ngknr)
  deallocate(igkignr)
  deallocate(occsvnr)
endif

return
end
