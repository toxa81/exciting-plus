subroutine wann_val(r,val)
use modmain
use mod_timer
implicit none
! arguments
real(8), intent(in) :: r(3)
complex(4), intent(out) :: val(nspinor,nwfplot)

integer ntr(3),n
integer is,ia,ias,ir0,io,lm,ig,i,j,l,ispn
real(8) ya(nprad),c(nprad),tp(2)
real(8) vr0(3),t1,r0,tr(3)
real(8) ur(0:lmaxvr,nrfmax)
complex(8) ylm(lmmaxvr)
complex(8) zt1(nspinor,nwfplot),zt2(nspinor,nwfplot),zt3
integer ikloc
real(8), external :: polynom 
logical, external :: vrinmt

zt1=zzero
zt2=zzero
if (vrinmt(r,is,ia,ntr,vr0,ir0,r0)) then
  call timer_start(1)
  ias=idxas(ia,is)
  call sphcrd(vr0,t1,tp)
  call genylm(lmaxvr,tp,ylm)
  tr(:)=ntr(1)*avec(:,1)+ntr(2)*avec(:,2)+ntr(3)*avec(:,3)
  do io=1,nrfmax
    do l=0,lmaxvr
      do j=1,nprad
        i=ir0+j-1
        ya(j)=urf(i,l,io,ias)
      end do
      ur(l,io)=polynom(0,nprad,spr(ir0,is),ya,c,r0)
    enddo !l
  enddo !io
  do ikloc=1,nkptloc
    zt3=exp(zi*dot_product(vkc(:,ikloc),tr(:)))
    do n=1,nwfplot
      do ispn=1,nspinor
        do io=1,nrfmax
          do lm=1,lmmaxvr
            zt1(ispn,n)=zt1(ispn,n)+zt3*wann_unkmt(lm,io,ias,ispn,n+firstwf-1,ikloc)* &
              ur(lm2l(lm),io)*ylm(lm)
          enddo
        enddo
      enddo
    enddo
  enddo
  !call zsync(zt1,nspinor*nwfplot,.true.,.false.)
  call timer_stop(1)
else
  call timer_start(2)
  do ikloc=1,nkptloc
    do ispn=1,nspinor
      do ig=1,ngk(1,ikloc)
        zt3=exp(zi*dot_product(vgkc(:,ig,1,ikloc),r(:)))/sqrt(omega)
        do n=1,nwfplot
          zt2(ispn,n)=zt2(ispn,n)+zt3*wann_unkit(ig,ispn,n+firstwf-1,ikloc)
        enddo
      end do
    enddo
  end do
  !call zsync(zt2,nspinor*nwfplot,.true.,.false.)
  call timer_stop(2)
endif
if (iproc.eq.0) then
  val(:,:)=(zt1(:,:)+zt2(:,:))/nkpt
endif
!call barrier(comm_world)
return
end
