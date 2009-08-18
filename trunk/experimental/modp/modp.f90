module modp

integer, parameter :: ntimers=100
real(8) :: timer(ntimers,2)

integer cart_ndims
integer, allocatable :: cart_dims(:)
integer, allocatable :: cart_x(:)
integer cart_ncomm
integer, allocatable :: cart_comm(:)

contains
!
! timer routines
!
subroutine timer_start(n,reset)
implicit none
integer, intent(in) :: n
logical, optional, intent(in) :: reset
real(8) cpu0
if (present(reset)) then
  if (reset) call timer_reset(n)
endif
call cpu_time(cpu0)
timer(n,1)=cpu0
return
end subroutine

subroutine timer_stop(n)
implicit none
integer, intent(in) :: n
real(8) cpu0
call cpu_time(cpu0)
timer(n,2)=timer(n,2)+cpu0-timer(n,1)
timer(n,1)=cpu0
return
end subroutine

subroutine timer_reset(n)
implicit none
integer, intent(in) :: n
timer(n,:)=0.d0
return
end subroutine

!
! mpi parallel routines
!
subroutine cart_initialize(cart_dims_)
use mpi
implicit none
! arguments
integer, dimension(:), intent(in) :: cart_dims_
! local variables
logical, allocatable :: cart_p(:)
integer ierr,i,j

cart_ndims=size(cart_dims_)
cart_ncomm=2**cart_ndims
if (allocated(cart_dims)) deallocate(cart_dims)
allocate(cart_dims(cart_ndims))
cart_dims=cart_dims_
if (allocated(cart_comm)) deallocate(cart_comm)
allocate(cart_comm(0:cart_ncomm))
cart_comm=MPI_COMM_NULL
if (allocated(cart_x)) deallocate(cart_x)
allocate(cart_x(cart_ndims))
cart_x=-1

allocate(cart_p(cart_ndims))
cart_p=.false.

call mpi_cart_create(MPI_COMM_WORLD,cart_ndims,cart_dims,cart_p, &
  .false.,cart_comm(0),ierr)
do i=1,cart_ncomm
  write(*,*)'configuration : ',i
  cart_p=.false.
  do j=1,cart_ndims
    if (iand(i-1,2**(j-1)).ne.0) cart_p(j)=.true.
  enddo
  call mpi_cart_sub(cart_comm(0),cart_p,cart_comm(i),ierr)
  write(*,*)cart_p,cart_comm(i)
enddo
return
end subroutine

end module