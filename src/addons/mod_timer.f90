module mod_timer

integer, parameter :: ntimers=100
real(8) :: timer_value(ntimers,2)
integer :: timer_count(ntimers)

contains

subroutine timestamp(fout,txt)
implicit none
integer, intent(in) :: fout
character(*), optional, intent(in) :: txt
integer values(8)
call date_and_time(values=values)
write(fout,*)
write(fout,'("timestamp: ",I2.2,".",I2.2,".",I4.4,"  ",I2.2,":",I2.2,":",I2.2)') &
  values(3),values(2),values(1),values(5),values(6),values(7)
if (present(txt)) write(fout,'("comment: ",A)')trim(adjustl(txt))
return
end subroutine


real(8) function cpu_seconds()
implicit none
integer values(8)
call date_and_time(values=values)
cpu_seconds=(values(3)*24+values(5))*3600.d0+values(6)*60.d0+&
  values(7)*1.d0+values(8)/1000.d0
return
end function

subroutine timer_start(n,reset)
implicit none
integer, intent(in) :: n
logical, optional, intent(in) :: reset
real(8) cpu0
if (present(reset)) then
  if (reset) call timer_reset(n)
endif
cpu0=cpu_seconds()
timer_value(n,1)=cpu0
return
end subroutine

subroutine timer_stop(n)
implicit none
integer, intent(in) :: n
real(8) cpu0
cpu0=cpu_seconds()
timer_value(n,2)=timer_value(n,2)+cpu0-timer_value(n,1)
timer_value(n,1)=cpu0
timer_count(n)=timer_count(n)+1
return
end subroutine

subroutine timer_reset(n)
implicit none
integer, intent(in) :: n
if (n.eq.0) then  
  timer_value(2:,:)=0.d0
  timer_count(2:)=0
else
  timer_value(n,:)=0.d0
  timer_count(n)=0
endif
return
end subroutine

real(8) function timer_get_value(n,average)
implicit none
integer, intent(in) :: n
logical, optional, intent(in) :: average
logical average_
integer i
average_=.false.
if (present(average)) average_=average
if (average_) then
  i=timer_get_count(n)
  if (i.ne.0) then
    timer_get_value=timer_value(n,2)/dble(i)
  else
    timer_get_value=-1.d0
  endif
else
  timer_get_value=timer_value(n,2)
endif
return
end function

integer function timer_get_count(n)
implicit none
integer, intent(in) :: n
timer_get_count=timer_count(n)
return
end function

end module

