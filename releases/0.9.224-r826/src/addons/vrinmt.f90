logical function vrinmt(vr,is,ia,tr,vr0,ir0,r0)
use modmain
implicit none
! arguments
real(8), intent(in) :: vr(3)
integer, intent(out) :: is
integer, intent(out) :: ia
integer, intent(out) :: tr(3)
real(8), intent(out) :: vr0(3)
integer, intent(out) :: ir0
real(8), intent(out) :: r0
! local variables
real(8) vr0l(3),r1(3),rmt2,pos(3)
integer i1,i2,i3,ir,np2

vrinmt=.false.
call getntr(avec,vr,tr,vr0l)
r1(:)=vr0l(1)*avec(:,1)+vr0l(2)*avec(:,2)+vr0l(3)*avec(:,3)

np2=nprad/2
do is=1,nspecies
  rmt2=rmt(is)**2
  do ia=1,natoms(is)
    do i1=-1,1
    do i2=-1,1
    do i3=-1,1
      pos(:)=atposc(:,ia,is)+i1*avec(:,1)+i2*avec(:,2)+i3*avec(:,3)
      vr0(:)=r1(:)-pos(:)
      r0=vr0(1)**2+vr0(2)**2+vr0(3)**2
      if (r0.lt.rmt2) then
        r0=sqrt(r0)
        do ir=1,nrmt(is)
          if (spr(ir,is).ge.r0) then
            if (ir.le.np2) then
              ir0=1
            else if (ir.gt.nrmt(is)-np2) then
              ir0=nrmt(is)-nprad+1
            else
              ir0=ir-np2
            end if
            r0=max(r0,spr(1,is))
            tr(1)=tr(1)+i1
            tr(2)=tr(2)+i2
            tr(3)=tr(3)+i3
            vrinmt=.true.
            return
          end if
        end do !ir
      end if
    end do
    end do
    end do
  end do
end do
return
end
