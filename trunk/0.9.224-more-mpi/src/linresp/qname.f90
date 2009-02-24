subroutine qname(ivq,name)
implicit none
integer, intent(in) :: ivq(3)
character*(*), intent(out) :: name
character*6 fmt1,str1
integer i
name='q_'
do i=1,3
  if (ivq(i).lt.0) then 
    fmt1='(I4.3)'
  else
    fmt1='(I3.3)'
  endif
  write(str1,fmt1)ivq(i)
  name=trim(name)//trim(adjustl(str1))
  if (i.lt.3) name=trim(name)//"_"
enddo
return
end

