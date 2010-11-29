subroutine checkherm(n,m,irow,jcol)
implicit none
integer, intent(in) :: n
complex(8), intent(in) :: m(n,n)
integer, intent(out) :: irow
integer, intent(out) :: jcol
real(8), parameter :: epsherm=1d-10
integer i,j
real(8) t1,t2
irow=-1
jcol=-1
t2=-1.d0
do i=1,n
  do j=1,n
    t1=abs(m(i,j)-dconjg(m(j,i)))
    if (t1.gt.epsherm.and.t1.gt.t2) then
      t2=t1
      irow=i
      jcol=j
    endif
  enddo
enddo
return
end