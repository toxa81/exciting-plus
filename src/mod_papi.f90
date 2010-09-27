module mod_papi
use mod_mpi_grid

integer, parameter :: papi_ntimers=100
integer :: papi_ncounters
integer :: papi_eventset
integer(8) :: papi_timer(2,papi_ntimers)
integer(8), allocatable :: papi_counter(:,:,:)
character*256, allocatable :: papi_events(:)
integer clockrate

contains

!-----------------!
! papi_initialize !
!-----------------!
subroutine papi_initialize(nevents,events)
implicit none
#ifdef _PAPI_
include 'f90papi.h'
#endif
! arguments
integer, intent(in) :: nevents
character*(*), intent(in) :: events(*)
#ifdef _PAPI_
! local variables
integer check,eventcode,i,ierr
character*256 errstr

check=PAPI_VER_CURRENT
call PAPIF_library_init(check)
if (check.ne.PAPI_VER_CURRENT) then
  write(*,'("Error: PAPI library version is out of date")')
  call pstop
endif
call PAPIF_get_clockrate(clockrate)
! get number of hardware counters
call PAPIF_num_counters(papi_ncounters)
! take minimum value
papi_ncounters=min(nevents,papi_ncounters)
allocate(papi_counter(papi_ncounters,2,papi_ntimers))
allocate(papi_events(papi_ncounters))
! create set of counters
papi_eventset=PAPI_NULL
call PAPIF_create_eventset(papi_eventset,check)
if (check.ne.PAPI_OK) then
  write(*,'("Error: PAPIF_create_eventset returned : ",I6)')check
  call pstop
endif
do i=1,papi_ncounters
  papi_events(i)=trim(adjustl(events(i)))
  call PAPIF_event_name_to_code(trim(adjustl(events(i))),eventcode,check)
  if (check.ne.PAPI_OK) then
    write(*,'("Error: PAPIF_event_name_to_code : ",I6)')check
    write(*,'("  event name : ",A)')trim(adjustl(events(i)))
    call pstop
  endif
  call PAPIF_add_event(papi_eventset,eventcode,check)  
  if (check.ne.PAPI_OK) then
    write(*,'("Error: PAPIF_add_event returned : ",I6)')check
    call PAPIF_perror(check,errstr,ierr)
    write(*,'(" Error string : : ",A)')trim(adjustl(errstr))   
    call pstop
  endif
enddo
papi_timer=0
papi_counter=0
if (papi_ncounters.eq.0) return
! start counters
call PAPIF_start(papi_eventset,check)
if (check.ne.PAPI_OK) then
  write(*,'("Error: PAPIF_start returned : ",I6)')check
  call pstop
endif
#endif
end subroutine

subroutine proc_load(a,b,c)
implicit none
real(8), intent(in) :: a
real(8), intent(in) :: b
real(8), intent(out) :: c
c=c+a*b
end subroutine

!---------------!
! papi_finalize !
!---------------!
subroutine papi_finalize
implicit none
#ifdef _PAPI_
include 'f90papi.h'
call PAPIF_shutdown 
#endif
end subroutine

!------------------!
! papi_timer_start !
!------------------!
subroutine papi_timer_start(n,reset)
implicit none
#ifdef _PAPI_
include 'f90papi.h'
#endif
! arguments
integer, intent(in) :: n
logical, optional, intent(in) :: reset
#ifdef _PAPI_
integer check
if (present(reset)) then
  if (reset) call papi_timer_reset(n)
endif
!call PAPIF_get_real_usec(papi_timer(1,n),check)
call PAPIF_get_real_cyc(papi_timer(1,n),check)
if (papi_ncounters.eq.0) return
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
implicit none
#ifdef _PAPI_
include 'f90papi.h'
#endif
! arguments
integer, intent(in) :: n
#ifdef _PAPI_
integer check,i
integer(8) counter0(papi_ncounters)
integer(8) time0
!call PAPIF_get_real_usec(time0,check)
call PAPIF_get_real_cyc(time0,check)
if (papi_ncounters.ne.0) then
  call PAPIF_read(papi_eventset,counter0,check)
  if (check.ne.PAPI_OK) then
    write(*,'("Error: PAPIF_read returned : ",I6)')check
    call pstop
  endif
endif
papi_timer(2,n)=papi_timer(2,n)+time0-papi_timer(1,n)
do i=1,papi_ncounters
  papi_counter(i,2,n)=papi_counter(i,2,n)+counter0(i)-papi_counter(i,1,n)
enddo
#endif
end subroutine

!------------------!
! papi_timer_reset !
!------------------!
subroutine papi_timer_reset(n)
implicit none
#ifdef _PAPI_
include 'f90papi.h'
#endif
! arguments
integer, intent(in) :: n
#ifdef _PAPI_
papi_timer(:,n)=0
papi_counter(:,:,n)=0
#endif
end subroutine

subroutine papi_report(fout,values,comment)
implicit none
integer, intent(in) :: fout
integer*8, intent(in) :: values(0:papi_ncounters)
character*(*), intent(in) :: comment
#ifdef _PAPI_
! local variables
integer i
integer(8) cycles
real(8) time
integer(8) fp_ops,l1_dca,l1_dcm

fp_ops=-1
l1_dca=-1
l1_dcm=-1
cycles=values(0)
time=cycles/(clockrate*1.0d6)
write(fout,'(60("-"))')
write(fout,'("PAPI report : ",A)')trim(adjustl(comment))
write(fout,'(60("-"))')
write(fout,'("approx. time (seconds) : ",G18.10)')time
write(fout,'("cycles : ",I20)')cycles
do i=1,papi_ncounters
  write(fout,'(A," : ",I20)')trim(adjustl(papi_events(i))),values(i)
  if (trim(adjustl(papi_events(i))).eq."PAPI_FP_OPS") fp_ops=values(i)
  if (trim(adjustl(papi_events(i))).eq."PAPI_L1_DCA") l1_dca=values(i)
  if (trim(adjustl(papi_events(i))).eq."PAPI_L1_DCM") l1_dcm=values(i)
enddo
write(fout,'(60("-"))')
write(fout,'("Derived metrics:")')
if (fp_ops.ge.0) then
  write(fout,'("Computational intensity (ops/cycle) : ",F8.4)')dble(fp_ops)/cycles
  write(fout,'("Performance (GFlops) : ",F8.4)')dble(fp_ops)/time/1.0d9  
endif
if (l1_dca.ge.0.and.l1_dcm.ge.0) then
  write(fout,'("D1 cache hit,miss ratio (%) : ",2F8.4)')100.d0*dble(l1_dca)/(l1_dca+l1_dcm),&
    100.d0*dble(l1_dcm)/(l1_dca+l1_dcm)
  if (l1_dcm.gt.0) then 
    write(fout,'("D1 cache utilization (hits/miss) : ",G18.10)')dble(l1_dca)/l1_dcm
  endif
endif
write(fout,'(60("-"))')
#endif
return
end subroutine

subroutine papi_timer_read(n,values)
implicit none
integer, intent(in) :: n
integer(8), intent(out) :: values(0:papi_ncounters)
#ifdef _PAPI_
integer i
values(0)=papi_timer(2,n)
do i=1,papi_ncounters
  values(i)=papi_counter(i,2,n)
enddo
#endif
return
end subroutine


end module

