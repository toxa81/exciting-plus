subroutine memcopy(src,dest,size)
implicit none
character, intent(in) :: src(*)
character, intent(out) :: dest(*)
integer, intent(in) :: size
dest(1:size)=src(1:size)
return
end