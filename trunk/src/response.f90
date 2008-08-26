subroutine response
use modmain
implicit none

integer, allocatable :: igishell(:)
integer, allocatable :: ishellng(:,:)

integer i,j,ngsh,ngshmin,ngvec1,ngvec_me,ngvec_chi
character*100 fname

real(8) vq0l(3)
integer vgq0l(3)

! initialise universal variables
call init0
call init1

!call writegshells

if (.not.spinpol) then
  spin_me=1
  spin_chi=1
endif

if (task.eq.400) fname='RESPONSE_ME.OUT'
if (task.eq.401) fname='RESPONSE_CHI0.OUT'
if (task.eq.402) fname='RESPONSE_CHI.OUT'

if (iproc.eq.0) then
  open(150,file=trim(fname),form='formatted',status='replace')
  if (ismpi) then
    write(150,'("Running in parallel mode on ",I4," proc.")')nproc
  else
    write(150,'("Running in serial mode")')
  endif
endif

if (task.eq.400) then
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
    write(*,'("Warning(response): minimum number of G-shells changed to: ",I4)')ngshmin
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
  
! calculate matrix elements
  do i=1,nvq0
    call response_me_v3(ivq0m_list(1,i),ngvec_me)
  enddo
  
  deallocate(igishell)
  deallocate(ishellng)
endif

if (task.eq.401) then
  do i=1,nvq0
    call response_chi0(ivq0m_list(1,i))
  enddo
endif

if (task.eq.402) then
  allocate(igishell(ngvec))
  allocate(ishellng(ngvec,2))
  call getgshells(ngsh,igishell,ishellng)
  
  ngvec_chi=ishellng(ngsh_chi,2)
! calculate chi
  do i=1,nvq0
    call response_chi(ivq0m_list(1,i),ngvec_chi)
  enddo
  
  deallocate(igishell)
  deallocate(ishellng)
endif

if (iproc.eq.0) close(150)

return
end
