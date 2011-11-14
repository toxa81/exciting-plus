subroutine writegshells
use modmain
implicit none

integer ngsh
integer ,allocatable :: igishell(:)
integer ,allocatable :: ishellng(:,:)

integer ig,ish

allocate(igishell(ngvec))
allocate(ishellng(ngvec,2))

call getgshells(ngsh,igishell,ishellng)
open(170,file='GSHELLS.OUT',form='formatted',status='replace')
write(170,*)
write(170,'("  G-vec.       lat.coord.      length(a.u.)   shell ")')
write(170,'(" ---------------------------------------------------")')
do ig=1,ishellng(ngsh,2)
  write(170,'(2X,I8,4X,3I5,4X,F12.6,5x,I4)')ig,ivg(:,ig),gc(ig),igishell(ig)
enddo
write(170,*)
write(170,'("  G-shell    num. G-vec in shell      total num. G-vec. ")')
write(170,'(" -------------------------------------------------------")')
do ish=1,ngsh
  write(170,'(2X,I4,14X,I4,18X,I8)')ish,ishellng(ish,1),ishellng(ish,2)
enddo
close(170)

deallocate(igishell,ishellng)

return
end
