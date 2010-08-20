module mod_papi
use mod_mpi_grid
!integer, parameter :: papi_nset=100
!integer(8) :: papi_eventset_values(10,papi_nset)
!real(8) :: papi_flops(papi_nset)

integer, parameter :: papi_ncounter=4
integer, parameter :: papi_ntimer=100
integer :: papi_eventset
integer(8) :: papi_timer(2,papi_ntimer)
integer(8) :: papi_counter(papi_ncounter,2,papi_ntimer)

contains

!-----------------!
! papi_initialize !
!-----------------!
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
call papi_start_set
papi_timer=0
papi_counter=0
#endif
papi_flops=-1000000.d0
end subroutine

!---------------!
! papi_finalize !
!---------------!
subroutine papi_finalize
#ifdef _PAPI_
implicit none
include 'f90papi.h'
!call papi_stop_set
call PAPIF_shutdown 
#endif
end subroutine

!----------------!
! papi_start_set !
!----------------!
subroutine papi_start_set
#ifdef _PAPI_
implicit none
include 'f90papi.h'
integer check,eventcode

papi_eventset=PAPI_NULL
call PAPIF_create_eventset(papi_eventset,check)
if (check.ne.PAPI_OK) then
  write(*,'("Error: PAPIF_create_eventset returned : ",I6)')check
  call pstop
endif

eventcode=PAPI_FP_OPS
call PAPIF_add_event(papi_eventset,eventcode,check)
if (check.ne.PAPI_OK) then
  write(*,'("Error: PAPIF_add_event returned : ",I6)')check
  call pstop
endif
!call PAPIF_get_real_usec(papi_eventset_values(1,iset),check)
call PAPIF_start(papi_eventset,check)
if (check.ne.PAPI_OK) then
  write(*,'("Error: PAPIF_start returned : ",I6)')check
  call pstop
endif
#endif
end subroutine

!------------------!
! papi_timer_start !
!------------------!
subroutine papi_timer_start(n)
#ifdef _PAPI_
implicit none
include 'f90papi.h'
integer, intent(in) :: n
logical, optional, intent(in) :: reset
integer check
if (present(reset)) then
  if (reset) call papi_timer_reset(n)
endif
call PAPIF_get_real_usec(papi_timer(1,n),check)
call PAPIF_read(papi_eventset,papi_counter(1,1,n),check)
if (check.ne.PAPI_OK) then
  write(*,'("Error: PAPIF_read returned : ",I6)')check
  call pstop
endif
#endif
end subroutine

!-----------------!
! papi_timer_stop !
!-----------------!
subroutine papi_timer_stop(n)
#ifdef _PAPI_
implicit none
include 'f90papi.h'
integer, intent(in) :: n
integer check,i
integer(8) counter0(papi_ncounter)
integer(8) time0
call PAPIF_get_real_usec(time0,check)
call PAPIF_read(papi_eventset,counter0,check)
if (check.ne.PAPI_OK) then
  write(*,'("Error: PAPIF_read returned : ",I6)')check
  call pstop
endif
papi_timer(2,n)=papi_timer(2,n)+time0-papi_timer(1,n)
do i=1,papi_ncounter
  papi_counter(2,i,n)=papi_counter(2,i,n)+counter0(i)-papi_counter(1,i,n)
enddo
#endif
end subroutine

!------------------!
! papi_timer_reset !
!------------------!
subroutine papi_timer_reset(n)
#ifdef _PAPI_
implicit none
include 'f90papi.h'
integer, intent(in) :: n
papi_timer(:,n)=0
papi_counter(:,:,n)=0
#endif
end subroutine

real(8) function papi_timer_get_value(i,n)
implicit none
integer, intent(in) :: i
integer, intent(in) :: n
real(8) t1,t2
t2=papi_timer(2,n)/1.0d6
t1=papi_counter(i,2,n)/t2
papi_timer_get_value=t1
end function

!---------------!
! papi_stop_set !
!---------------!
!subroutine papi_stop_set
!#ifdef _PAPI_
!implicit none
!include 'f90papi.h'
!integer check
!call PAPIF_stop(papi_eventset,papi_eventset_values(3,iset),check)
!if (check.ne.PAPI_OK) then
!  write(*,'("Error: PAPIF_stop returned : ",I6)')check
!  call pstop
!endif
!call PAPIF_get_real_usec(papi_eventset_values(2,iset),check)
!tottime=(papi_eventset_values(2,iset)-papi_eventset_values(1,iset))/1.0d6
!papi_flops(iset)=papi_eventset_values(3,iset)/tottime
!#endif
!end subroutine

!real(8) function papi_mflops(iset)
!implicit none
!integer, intent(in) :: iset
!papi_mflops=papi_flops(iset)/1.0d6
!end function
!

end module
