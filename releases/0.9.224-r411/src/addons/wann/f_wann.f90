subroutine f_wann(n,ispn,d1,itr,vrc,val)
use modmain
implicit none
integer, intent(in) :: n
integer, intent(in) :: ispn
real(8), intent(in) :: d1
integer, intent(in) :: itr(3)
real(8), intent(in) :: vrc(3)
real(8), intent(out) :: val(2)

integer ik,i1,i2,i3,np2,is,ia,ias,ir,ir0,io,l,j,m,i,lm,ig
integer ntr1(3)
real(8) t1,rmt2
real(8) dtr(3),trvec(3),vrc1(3),v1(3),vr0l(3),tp(2),pos(3)
complex(8) zt1,zt2,zt3,ylm(lmmaxvr)
logical l1
real(8) ya(nprad),c(nprad),t2
real(8), external :: polynom

dtr(:)=1.d0*itr(:)
call r3mv(avec,dtr,trvec)
vrc1(:)=vrc(:)-trvec(:)

t1=sqrt(vrc1(1)**2+vrc1(2)**2+vrc1(3)**2)
if (t1.gt.d1) then
  val(1)=0.d0
  val(2)=0.d0
  return
endif

np2=nprad/2
zt1=dcmplx(0.d0,0.d0)
zt2=dcmplx(0.d0,0.d0)

call getntr(avec,vrc1,ntr1,vr0l)

l1=.false.
! check if point is in a muffin-tin
do is=1,nspecies
  rmt2=rmt(is)**2
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do i1=ntr1(1)-1,ntr1(1)+1
    do i2=ntr1(2)-1,ntr1(2)+1
    do i3=ntr1(3)-1,ntr1(3)+1
      trvec(:)=i1*avec(:,1)+i2*avec(:,2)+i3*avec(:,3)
      pos(:)=atposc(:,ia,is)+trvec(:)
      v1(:)=vrc1(:)-pos(:)
      t1=v1(1)**2+v1(2)**2+v1(3)**2
      if (t1.lt.rmt2) then
        call sphcrd(v1,t1,tp)
        call genylm(lmaxvr,tp,ylm)
        do ir=1,nrmt(is)
          if (spr(ir,is).ge.t1) then
            if (ir.le.np2) then
              ir0=1
            else if (ir.gt.nrmt(is)-np2) then
              ir0=nrmt(is)-nprad+1
            else
              ir0=ir-np2
            end if
            t1=max(t1,spr(1,is))
            do io=1,nrfmax
              do l=0,lmaxvr
                do j=1,nprad
                  i=ir0+j-1
                  ya(j)=urf(i,l,io,ias)
                end do
                t2=polynom(0,nprad,spr(ir0,is),ya,c,t1)
                do m=-l,l
                  lm=idxlm(l,m)
                  zt3=dcmplx(0.d0,0.d0)
                  do ik=1,nkpt
                    zt3=zt3+wann_unkmt(lm,io,ias,n,ispn,ik)* &
                      exp(zi*dot_product(vkc(:,ik),trvec(:)))
                  enddo
                  zt1=zt1+zt3*t2*ylm(lm)
                end do !m
              end do  !l
            enddo !io
            l1=.true.
            goto 10
          end if
        end do !ir
      end if
    end do
    end do
    end do
  end do
end do !is
10 continue
! otherwise use interstitial function
if (.not.l1) then
  do ik=1,nkpt
    do ig=1,ngk(1,ik)
      t1=vgkc(1,ig,1,ik)*vrc1(1)+vgkc(2,ig,1,ik)*vrc1(2)+vgkc(3,ig,1,ik)*vrc1(3)
      zt2=zt2+exp(zi*t1)*wann_unkit(ig,n,ispn,ik)/sqrt(omega)
    enddo
  enddo
endif
val(1)=dreal(zt1+zt2)/nkpt
val(2)=dimag(zt1+zt2)/nkpt

return
end
