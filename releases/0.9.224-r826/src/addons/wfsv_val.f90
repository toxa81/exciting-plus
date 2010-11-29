subroutine wfsv_val(vrc,vpc,vgpc,ngp,wfsvmt,wfsvit,val)
use modmain
implicit none
real(8), intent(in) :: vrc(3)
real(8), intent(in) :: vpc(3)
real(8), intent(in) :: vgpc(3,ngkmax)
integer, intent(in) :: ngp
complex(8), intent(in) :: wfsvmt(lmmaxvr,nrfmax,natmtot,nspinor)
complex(8), intent(in) :: wfsvit(ngkmax,nspinor)
complex(8), intent(out) :: val(nspinor) 
real(8), external :: polynom 
logical, external :: vrinmt
integer is,ia,ias,ntr(3),ir0,i,j,l,lm,io,ispn,ig
real(8) vr0(3),r0,t1,tp(2),ya(nprad),c(nprad),ur(0:lmaxvr,nrfmax),tr(3)
complex(8) ylm(lmmaxvr),zt1

val=zzero
if (vrinmt(vrc,is,ia,ntr,vr0,ir0,r0)) then
  tr(:)=ntr(1)*avec(:,1)+ntr(2)*avec(:,2)+ntr(3)*avec(:,3)
  ias=idxas(ia,is)
  call sphcrd(vr0,t1,tp)
  call genylm(lmaxvr,tp,ylm)
  do io=1,nrfmax
    do l=0,lmaxvr
      do j=1,nprad
        i=ir0+j-1
        ya(j)=urf(i,l,io,ias)
      end do
      ur(l,io)=polynom(0,nprad,spr(ir0,is),ya,c,r0)
    enddo !l
  enddo !io
  zt1=exp(zi*dot_product(vpc(:),tr(:)))
  do ispn=1,nspinor
    do io=1,nrfmax
      do lm=1,lmmaxvr
        val(ispn)=val(ispn)+zt1*wfsvmt(lm,io,ias,ispn)*ur(lm2l(lm),io)*ylm(lm)
      enddo
    enddo
  enddo
else
  do ispn=1,nspinor
    do ig=1,ngp
      zt1=exp(zi*dot_product(vgpc(:,ig),vrc(:)))/sqrt(omega)
      val(ispn)=val(ispn)+zt1*wfsvit(ig,ispn)
    enddo
  enddo
endif
return
end