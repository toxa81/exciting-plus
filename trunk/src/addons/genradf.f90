subroutine genradf
use modmain
implicit none
! find the new linearisation energies
call timer_start(t_lin_en)
call linengy
call timer_stop(t_lin_en)
! generate the APW radial functions
call timer_start(t_apw_rad)
call genapwfr
! generate the local-orbital radial functions
call genlofr
call timer_stop(t_apw_rad)
! collect radial-muffint tin functions
call getufr
! compute the product of radial functions
call genufrp
return
end subroutine
