subroutine getntr(avec,vrcnr,ntr,vrl)
!use modmain
implicit none
real(8), intent(in) :: avec(3,3)
real(8), intent(in) :: vrcnr(3)
integer, intent(out) :: ntr(3)
real(8), intent(out) :: vrl(3)
real(8) a(3,3)
real(8) b(3)
integer i,j,ipiv(3)
real(8) work(200)
integer lwork
! r=r0+T=r0+n1*a1+n2*a2+n3*a3=f1*a1+f2*a2+f3*a3
! system of linear equtions for f1,f2,f3
! a1*r=a1*(f1*a1+f2*a2+f3*a3)
! a2*r=a2*(f1*a1+f2*a2+f3*a3)
! a3*r=a3*(f1*a1+f2*a2+f3*a3)
do j=1,3
  do i=1,j
    a(i,j)=dot_product(avec(:,i),avec(:,j))
  enddo
  b(j)=dot_product(avec(:,j),vrcnr(:))
enddo
lwork=200
call dsysv('U',3,1,a,3,ipiv,b,3,work,lwork,i)
do i=1,3
  ntr(i)=floor(b(i))
  vrl(i)=b(i)-ntr(i)
enddo
return
end
