subroutine getntr(vrcnr,ntr,vrl)
use modmain
implicit none
real(8), intent(in) :: vrcnr(3)
integer, intent(out) :: ntr(3)
real(8), intent(out) :: vrl(3)
real(8) f(3)
integer i
logical l1
! r=r0+T=r0+n1*a_1+n2*a_2+n3*a_3=f1*a_1+f2*a_2+f3*a_3
! f=a^{-1}*r
call r3mv(ainv,vrcnr,f)
l1=.false.
do i=1,3
  ntr(i)=floor(f(i))
  vrl(i)=f(i)-dble(ntr(i))
  if (vrl(i).eq.1.d0) then
    vrl(i)=0.d0
    ntr(i)=ntr(i)+1
  endif
  if (vrl(i).lt.0.d0.or.vrl(i).ge.1.d0) l1=.true.
enddo
if (l1) then
  write(*,*)
  write(*,'("Warning(getntr): wrong lattice coordinates for reduced&
    & cartesian vector")')
  write(*,'("  original vector : ",3G18.10)')vrcnr
  write(*,'("  lattice translations : ",3I4)')ntr
  write(*,'("  translation vector : ",3G18.10)')&
    avec(:,1)*ntr(1)+avec(:,2)*ntr(2)+avec(:,3)*ntr(3)  
  write(*,'("  lattice coordinates of reduced vector : ",3G18.10)')vrl
endif
return
end
