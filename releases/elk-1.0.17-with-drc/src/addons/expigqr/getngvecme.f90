subroutine getngvecme
use modmain
use mod_addons_q
implicit none
integer ngsh,gshq0,i,j
integer, allocatable :: igishell(:)
integer, allocatable :: ishellng(:,:)
allocate(igishell(ngvec))
allocate(ishellng(ngvec,2))
call getgshells(ngsh,igishell,ishellng)
ngvecme=ishellng(gqsh,2)
deallocate(igishell)
deallocate(ishellng)
return
end