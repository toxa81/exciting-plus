subroutine sic_gndstate
use modmain
implicit none
integer iter

!lmpi_grid=.true.
!mpi_grid=(/1,1,1/)
do iter=1,1
  if (iproc.eq.0) write(*,*)'SIC iteration ',iter
  call gndstate
  call genvsic
enddo

return
end