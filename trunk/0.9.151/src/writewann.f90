subroutine writewann
use modmain
use modwann
implicit none
integer ik,ispn,n,i

open(50,file='WFC.OUT',status='replace',form='formatted')
do ik=1,nkpt
  write(50,'("  ik : ",I3)')ik
  do n=1,wf_dim
    write(50,'("WF : ",I3)')n
    do i=1,nstfv
      write(50,'("    istfv : ",I3,"   a : ",2F12.6)')i,(abs(wfc(n,i,ispn,ik)),ispn=1,wann_nspins)
    enddo
  enddo
enddo
close(50)

return
end
