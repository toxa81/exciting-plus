subroutine sic_localize
use modmain
use mod_sic
implicit none
integer i

sic=.true.
do i=1,15
  call sic_main
  call sic_update_umtrx
enddo
return
end subroutine

