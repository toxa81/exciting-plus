subroutine qname(ivq,name)
implicit none
integer, intent(in) :: ivq(3)
character*(*), intent(out) :: name
character*6 str1
integer i
name='q_'
do i=1,3
  write(str1,'(I6)')ivq(i)
  name=trim(name)//trim(adjustl(str1))
  if (i.lt.3) name=trim(name)//"_"
enddo
return
end

