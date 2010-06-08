subroutine sic_gndstate
use modmain
implicit none
integer iter

do iter=1,2
  if (iproc.eq.0) write(*,*)'SIC iteration ',iter
  call gndstate
  call mpi_world_barrier
  call genvsic
  call mpi_world_barrier
enddo

return
end