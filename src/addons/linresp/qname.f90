subroutine getqname(ivq,name)
implicit none
integer, intent(in) :: ivq(3)
character*(*), intent(out) :: name
character*6 str1
integer i
name="q_"
do i=1,3
  write(str1,'(I6)')ivq(i)
  name=trim(name)//trim(adjustl(str1))
  if (i.lt.3) name=trim(name)//"_"
enddo
return
end

subroutine getqdir(iq,ivq,name)
implicit none
integer, intent(in) :: iq
integer, intent(in) :: ivq(3)
character*(*), intent(out) :: name
character*6 str1
call getqname(ivq,name)
write(str1,'(I4.4)')iq
name="./q/iq_"//trim(adjustl(str1))//"__"//trim(name)
return
end

