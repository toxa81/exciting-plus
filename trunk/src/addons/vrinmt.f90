logical function vrinmt(vrcnr,is,ia,ntr,vrc0,ir0,r0)
use modmain
implicit none
! arguments
real(8), intent(in) :: vrcnr(3)
integer, intent(out) :: is
integer, intent(out) :: ia
integer, intent(out) :: ntr(3)
real(8), intent(out) :: vrc0(3)
integer, intent(out) :: ir0
real(8), intent(out) :: r0
! local variables
real(8) vrl(3),vrc(3),rmt2,pos(3)
integer i1,i2,i3,ir,np2

vrinmt=.false.
call getntr(vrcnr,ntr,vrl)
call r3mv(avec,vrl,vrc)

np2=nprad/2
do is=1,nspecies
  rmt2=rmt(is)**2
  do ia=1,natoms(is)
    do i1=-1,1
      do i2=-1,1
        do i3=-1,1
          pos(:)=atposc(:,ia,is)+i1*avec(:,1)+i2*avec(:,2)+i3*avec(:,3)
          vrc0(:)=vrc(:)-pos(:)
          r0=vrc0(1)**2+vrc0(2)**2+vrc0(3)**2
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
                ntr(1)=ntr(1)+i1
                ntr(2)=ntr(2)+i2
                ntr(3)=ntr(3)+i3
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
