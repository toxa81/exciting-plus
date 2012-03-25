module mod_util

contains

real(8) function rintegrate(nr,r,f,m,g)
implicit none
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
real(8), intent(in) :: f(nr)
integer, optional, intent(in) :: m
real(8), optional, intent(out) :: g(nr)
!
real(8), allocatable :: f0(:),g0(:),cf(:,:)
integer i,m0
!
! r^{m} weight
m0=2
if (present(m)) m0=m
allocate(f0(nr),g0(nr),cf(4,nr))
do i=1,nr
  f0(i)=f(i)*(r(i)**m0)
enddo
call fderiv(-1,nr,r,f0,g0,cf)
rintegrate=g0(nr)
if (present(g)) g(:)=g0(:)
deallocate(f0,g0,cf)
return
end function

complex(8) function zintegrate(nr,r,f,m,g)
implicit none
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
complex(8), intent(in) :: f(nr)
integer, optional, intent(in) :: m
complex(8), optional, intent(out) :: g(nr)
!
real(8), allocatable :: f0r(:),f0i(:),g0r(:),g0i(:),cf(:,:)
integer i,m0
!
! r^{m} weight
m0=2
if (present(m)) m0=m
allocate(f0r(nr),f0i(nr),g0r(nr),g0i(nr),cf(4,nr))
do i=1,nr
  f0r(i)=dreal(f(i))*(r(i)**m0)
  f0i(i)=dimag(f(i))*(r(i)**m0)
enddo
call fderiv(-1,nr,r,f0r,g0r,cf)
call fderiv(-1,nr,r,f0i,g0i,cf)
zintegrate=dcmplx(g0r(nr),g0i(nr))
if (present(g)) g(:)=dcmplx(g0r(:),g0i(:))
deallocate(f0r,f0i,g0r,g0i,cf)
return
end function

end module
