! DJBHash function
! http://www.partow.net/programming/hashfunctions/index.html
integer function hash(str,len)
implicit none
character, intent(in) :: str(*)
integer, intent(in) :: len
integer i
hash=5381
do i=1,len
  hash=hash*32+hash+ichar(str(i))
enddo
end function