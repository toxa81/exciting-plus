subroutine sic_test_vme
use modmain
use mod_sic
implicit none
!
sic=.true.
! initialise universal variables
call init0
call init1
if (.not.mpi_grid_in()) return
wproc=mpi_grid_root()

call sic_read_data(.false.)

call sic_genvme_dotp(.false.)

return
end subroutine
