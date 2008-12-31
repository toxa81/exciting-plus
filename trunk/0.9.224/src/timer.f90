subroutine timer_start(n)
use modmain
implicit none
integer, intent(in) :: n
real(8) cpu0

call cpu_time(cpu0)
timer(n,1)=cpu0
return
end

subroutine timer_stop(n)
use modmain
implicit none
integer, intent(in) :: n
real(8) cpu0

call cpu_time(cpu0)
timer(n,2)=timer(n,2)+cpu0-timer(n,1)
timer(n,1)=cpu0
return
end

subroutine timer_reset(n)
use modmain
implicit none
integer, intent(in) :: n
timer(n,:)=0.d0
return
end