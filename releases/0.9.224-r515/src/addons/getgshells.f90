subroutine getgshells(ngsh,igishell,ishellng)
use modmain
implicit none
! arguments
integer, intent(out) :: ngsh
integer, intent(out) :: igishell(ngvec)
integer, intent(out) :: ishellng(ngvec,2)

integer ish,ig

ish=1
ig=0
do while (ig.lt.ngvec)
  ig=ig+1
  igishell(ig)=ish
  if (abs(gc(ig+1)-gc(ig)).gt.epslat) then
    ish=ish+1
  endif
enddo 

ngsh=ish-1

do ish=1,ngsh
  ishellng(ish,1)=0
  do ig=1,ngvec
    if (igishell(ig).eq.ish) ishellng(ish,1)=ishellng(ish,1)+1
  enddo
  ishellng(ish,2)=sum(ishellng(1:ish,1))
enddo
 
return
end

