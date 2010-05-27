real(8) function rfinteg(nr,r,f)
implicit none
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
real(8), intent(in) :: f(nr)
real(8) t1
integer ir
t1=0.d0
do ir=1,nr-1
  t1=t1+0.5d0*(f(ir)+f(ir+1))*(r(ir+1)-r(ir))
enddo
rfinteg=t1
return
end

