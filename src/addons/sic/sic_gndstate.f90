subroutine sic_gndstate
use modmain
implicit none

do iitersic=1,nitersic
  if (iproc.eq.0) write(*,*)'SIC iteration ',iitersic
  !call gndstate
  !call mpi_world_barrier
  call sic_genvwan
  call mpi_world_barrier
enddo

return
end