program writeinput
use modmain
implicit none
logical lexist

inquire(file='exciting.in',exist=lexist)
if (lexist) then
  call readinput
else
  avec=0.d0
  avec(1,1)=1.d0
  avec(2,2)=1.d0
  avec(3,3)=1.d0
endif

open(50,file='exciting.out',form='formatted',status='replace')
write(50,'("tasks")')
write(50,'(2x,I1)')1
write(50,*)
write(50,'("avec")')
write(50,'(3G18.10)')avec(:,1)
write(50,'(3G18.10)')avec(:,2)
write(50,'(3G18.10)')avec(:,3)
close(50)






end
