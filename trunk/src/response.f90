subroutine response
use modmain
implicit none

! initialise universal variables
call init0
call init1

if (iproc.eq.0) then
  open(150,file='RESPONSE.OUT',form='formatted',status='replace')
  if (ismpi) then
    write(150,'("Running in parallel mode on ",I4," proc.")')nproc
  else
    write(150,'("Running in serial mode")')
  endif    
endif

if (task.eq.400) then
  call response_me
endif

if (task.eq.401) then
  call response_chi0
endif

if (task.eq.402) then
  call response_chi
endif

if (iproc.eq.0) close(150)

return
end


subroutine response_chi
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none

return
end


