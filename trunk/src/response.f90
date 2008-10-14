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

integer i,j,ngsh,ngshmin,ngvec1,ngvec_me,ngvec_chi,ik,ikloc,ig, &
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
  ngshmin=1
  do i=1,nvq0
! q-vector in lattice coordinates
    do j=1,3
      vq0l(j)=1.d0*ivq0m_list(j,i)/ngridk(j)
    enddo
! find G-vector which brings q0 to first BZ
    vgq0l(:)=floor(vq0l(:))
    ngshmin=max(ngshmin,igishell(ivgig(vgq0l(1),vgq0l(2),vgq0l(3))))
  enddo !i
  if (ngshmin.gt.ngsh_me) then
    write(*,*)
    write(*,'("Warning(response): number of G-shells was changed from ",&
      &I4," to ",I4)')ngsh_me,ngshmin
    write(*,*)
    ngsh_me=ngshmin
  endif
! test if G-shell is closed
  ngvec_me=ishellng(ngsh_me,2)
  if (abs(gc(ngvec_me)-gc(ngvec_me+1)).lt.epslat) then
    write(*,*)
    write(*,'("Bug(response): G-shell is not closed")')
    write(*,*)
    call pstop
  endif
  deallocate(igishell)
  deallocate(ishellng)
endif

if (task.eq.402.or.task.eq.403) then
  allocate(igishell(ngvec))
  allocate(ishellng(ngvec,2))
  call getgshells(ngsh,igishell,ishellng)
  ngvec_chi=ishellng(ngsh_chi,2)
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
    write(150,'("Size of acoeff array (Mb) : ",I6)')lmmaxvr*nrfmax*natmtot*nstsv*nkptlocnr(0)/1024/1024
    write(150,'("Size of evec array (Mb)   : ",I6)')nmatmax*nstsv*nkptlocnr(0)/1024/1024
    write(150,*)
    write(150,'("Reading eigen-vectors")')
    call flushifc(150)
  endif
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
    endif
  enddo !ikloc
  if (iproc.eq.0) then
    write(150,'("Done.")')
    call flushifc(150)
  endif
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
    call response_me(ivq0m_list(1,i),ngvec_me,ikptlocnr,nkptlocnr, &
      wfsvmtloc,wfsvitloc,ngknr,igkignr,occsvnr)
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
    call response_chi(ivq0m_list(1,i),ngvec_chi)
  enddo
endif

if (task.eq.403) then
  do i=1,nvq0
    call response_me(ivq0m_list(1,i),ngvec_me,ikptlocnr,nkptlocnr, &
      wfsvmtloc,wfsvitloc,ngknr,igkignr,occsvnr)
    call response_chi0(ivq0m_list(1,i),ikptlocnr,nkptlocnr)
    call response_chi(ivq0m_list(1,i),ngvec_chi)
  enddo
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
