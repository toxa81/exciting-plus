subroutine response
use modmain
implicit none

integer i,j,ngsh,ngshmax,ngvec1

real(8) vq0l(3)
integer vgq0l(3)

! initialise universal variables
call init0
call init1

if (iproc.eq.0) then
  open(150,file='RESPONSE.OUT',form='formatted',status='replace')
  if (ismpi) then
    write(150,'("Running in parallel mode on ",I4," proc.")')nproc
  else
    write(150,'("Running in serial mode")')
  endif
endif

if (task.eq.400) then
! find minimum number of G-shells
  ngshmax=1
  do i=1,nvq0
! q-vector in lattice coordinates
    do j=1,3
      vq0l(j)=1.d0*ivq0m_list(j,i)/ngridk(j)
    enddo
! find G-vector which brings q0 to first BZ
    vgq0l(:)=floor(vq0l(:))
    ngsh=1
 10 continue
    call getngvec(ngsh,ngvec1)
    if (ivgig(vgq0l(1),vgq0l(2),vgq0l(3)).gt.ngvec1) then
      ngsh=ngsh+1
      goto 10
    endif
    ngshmax=max(ngsh,ngshmax)
  enddo !i
  if (ngshmax.gt.ngsh_me) then
    write(*,*)
    write(*,'("Error(response): minimum number of G-shells must be: ",I4)')ngshmax
    write(*,*)
    call pstop
  endif
! calculate matrix elements
  do i=1,nvq0
    call response_me(ivq0m_list(1,i))
  enddo
endif

if (task.eq.401) then
  do i=1,nvq0
    call response_chi0(ivq0m_list(1,i))
  enddo
endif

if (task.eq.402) then
  do i=1,nvq0
    call response_chi(ivq0m_list(1,i))
  enddo
endif

if (iproc.eq.0) close(150)

return
end
