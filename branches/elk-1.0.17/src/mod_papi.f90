module mod_papi
use mod_mpi_grid
integer, parameter :: papi_nset=100
integer :: papi_eventset(papi_nset)
integer(8) :: papi_eventset_values(10,papi_nset)
real(8) :: papi_flops(papi_nset)
contains

subroutine papi_initialize
#ifdef _PAPI_
implicit none
include 'f90papi.h'
integer check
check=PAPI_VER_CURRENT
call PAPIF_library_init(check)
if (check.ne.PAPI_VER_CURRENT) then
  write(*,'("Error: PAPI library version is out of date")')
  call pstop
endif
#endif
papi_flops=-1000000.d0
end subroutine

subroutine papi_start_set(iset)
#ifdef _PAPI_
implicit none
include 'f90papi.h'
integer, intent(in) :: iset
integer check,eventcode

papi_eventset(iset)=PAPI_NULL
call PAPIF_create_eventset(papi_eventset(iset),check)
if (check.ne.PAPI_OK) then
  write(*,'("Error: PAPIF_create_eventset returned : ",I6)')check
  call pstop
endif

eventcode=PAPI_FP_OPS
call PAPIF_add_event(papi_eventset(iset),eventcode,check)
if (check.ne.PAPI_OK) then
  write(*,'("Error: PAPIF_add_event returned : ",I6)')check
  call pstop
endif

call PAPIF_get_real_usec(papi_eventset_values(1,iset),check)
call PAPIF_start(papi_eventset(iset),check)
if (check.ne.PAPI_OK) then
  write(*,'("Error: PAPIF_start returned : ",I6)')check
  call pstop
endif
#endif
end subroutine

subroutine papi_stop_set(iset)
#ifdef _PAPI_
implicit none
include 'f90papi.h'
integer, intent(in) :: iset
integer check
real(8) tottime
call PAPIF_stop(papi_eventset(iset),papi_eventset_values(3,iset),check)
if (check.ne.PAPI_OK) then
  write(*,'("Error: PAPIF_stop returned : ",I6)')check
  call pstop
endif
call PAPIF_get_real_usec(papi_eventset_values(2,iset),check)
tottime=(papi_eventset_values(2,iset)-papi_eventset_values(1,iset))/1.0d6
papi_flops(iset)=papi_eventset_values(3,iset)/tottime
#endif
end subroutine

real(8) function papi_mflops(iset)
implicit none
integer, intent(in) :: iset
papi_mflops=papi_flops(iset)/1.0d6
end function


end module
