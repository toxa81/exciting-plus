subroutine f_wann(n,ispn,d1,itr,vrc,val)
use modmain
use mod_wannier
implicit none
integer, intent(in) :: n
integer, intent(in) :: ispn
real(8), intent(in) :: d1
integer, intent(in) :: itr(3)
real(8), intent(in) :: vrc(3)
real(8), intent(out) :: val(2)

integer is,ia,ias,ir0,io,l,j,i,lm,ig
integer ntr(3),ikloc
real(8) tr(3),trvec(3),vrc1(3),tp(2),vr0(3),r0,t1
real(8) ur(0:lmaxvr,nufrmax)
complex(8) zt1,zt2,zt3,ylm(lmmaxvr)
real(8) ya(nprad),c(nprad)
real(8), external :: polynom
logical, external :: vrinmt

tr(:)=1.d0*itr(:)
call r3mv(avec,tr,trvec)
vrc1(:)=vrc(:)-trvec(:)

r0=sqrt(vrc1(1)**2+vrc1(2)**2+vrc1(3)**2)
if (r0.gt.d1) then
  val(1)=0.d0
  val(2)=0.d0
  return
endif

zt1=zzero
zt2=zzero
if (vrinmt(vrc1,is,ia,ntr,vr0,ir0,r0)) then
  ias=idxas(ia,is)
  call sphcrd(vr0,t1,tp)
  call genylm(lmaxvr,tp,ylm)
  tr(:)=ntr(1)*avec(:,1)+ntr(2)*avec(:,2)+ntr(3)*avec(:,3)
  do io=1,nufrmax
    do l=0,lmaxvr
      do j=1,nprad
        i=ir0+j-1
        ya(j)=ufr(i,l,io,ias2ic(ias))
      end do
      ur(l,io)=polynom(0,nprad,spr(ir0,is),ya,c,r0)
    enddo !l
  enddo !io
  do ikloc=1,nkpt
    zt3=exp(zi*dot_product(vkc(:,ikloc),tr(:)))
    do io=1,nufrmax
      do lm=1,lmmaxvr
        zt1=zt1+zt3*wann_unkmt(lm,io,ias,ispn,n,ikloc)* &
          ur(lm2l(lm),io)*ylm(lm)
      enddo
    enddo
  enddo
else
  do ikloc=1,nkpt
    do ig=1,ngk(1,ikloc)
      zt3=exp(zi*dot_product(vgkc(:,ig,1,ikloc),vrc1(:)))/sqrt(omega)
      zt2=zt2+zt3*wann_unkit(ig,ispn,n,ikloc)
    enddo
  enddo
endif
val(1)=dreal(zt1+zt2)/nkpt
val(2)=dimag(zt1+zt2)/nkpt

return
end
