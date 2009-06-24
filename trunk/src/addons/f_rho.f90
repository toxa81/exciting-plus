subroutine f_rho(vrc,val)
use modmain
implicit none
real(8), intent(in) :: vrc(3)
real(8), intent(out) :: val

integer ntr(3)
real(8) vr0l(3)

call getntr(avec,vrc,ntr,vr0l)
call rfarray(lmaxvr,lmmaxvr,rhomt,rhoir,1,vr0l,val)

return
end
