! DJBHash function
! http://www.partow.net/programming/hashfunctions/index.html
integer function hash(str,n)
implicit none
integer, intent(in) :: n
character, intent(in) :: str(n)
integer i
hash=5381
do i=1,n
  hash=hash*32+hash+ichar(str(i))
enddo
end function