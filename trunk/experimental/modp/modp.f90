module modp

integer, parameter :: ntimers=100
real(8) :: timer(ntimers,2)

integer nproc
data nproc/1/
integer iproc
data iproc/0/
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
subroutine mpi_world_initialize
use mpi
implicit none
integer ierr
call mpi_init(ierr)
call mpi_comm_size(MPI_COMM_WORLD,nproc,ierr)
call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierr)
return
end subroutine

subroutine mpi_world_finalize
use mpi
implicit none
integer ierr
call mpi_barrier(MPI_COMM_WORLD,ierr)
call mpi_finalize(ierr)
return
end subroutine


subroutine mpi_cart_initialize(cart_dims_)
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

if (cart_ingrid()) then
  call mpi_cart_get(cart_comm(0),cart_ndims,cart_dims,cart_p,cart_x,ierr)
  do i=1,cart_ncomm
    cart_p=.false.
    do j=1,cart_ndims
      if (iand(i-1,2**(j-1)).ne.0) cart_p(j)=.true.
    enddo
    call mpi_cart_sub(cart_comm(0),cart_p,cart_comm(i),ierr)
  enddo
!  if (cart_root()) then
!    write(*,'("MPI grid dimensions : ",10I4)')cart_dims
!  endif
endif !cart_ingrid()
return
end subroutine

subroutine mpi_cart_finalize
use mpi
implicit none
integer i,ierr

if (cart_ingrid()) then
  do i=0,cart_ncomm
    call mpi_comm_free(cart_comm(i),ierr)
  enddo
endif
deallocate(cart_dims)
deallocate(cart_x)
deallocate(cart_comm)
return
end subroutine

logical function cart_ingrid()
use mpi
implicit none
cart_ingrid=(cart_comm(0).ne.MPI_COMM_NULL)
return
end function

logical function cart_root(dims)
implicit none
! arguments
integer, optional, dimension(:), intent(in) :: dims
! local variables
integer i
logical l1
if (.not.present(dims)) then
  if (all(cart_x.eq.0)) then
    cart_root=.true.
  else
    cart_root=.false.
  endif
else
  l1=.true.
  do i=1,size(dims)
    if (cart_x(dims(i)).ne.0) l1=.false.
  enddo
  cart_root=l1
endif
return
end function

logical function cart_side(dims)
implicit none
! arguments
integer, dimension(:), intent(in) :: dims
! local variables
integer dims1(cart_ndims)
integer i
logical l1
dims1=1
do i=1,size(dims)
  dims1(dims(i))=0
enddo
l1=.true.
do i=1,cart_ndims
  if (dims1(i).eq.1.and.cart_x(i).ne.0) l1=.false.
enddo
cart_side=l1
return
end function

integer function get_cart_comm(dims)
implicit none
! arguments
integer, optional, dimension(:), intent(in) :: dims
! local variables
integer i,j

if (.not.present(dims)) then
  get_cart_comm=cart_comm(0)
else
  j=1
  do i=1,size(dims)
    j=j+2**(dims(i)-1)
  enddo
  get_cart_comm=cart_comm(j)
endif
return
end function

subroutine d_cart_bcast(val,n,dims,side)
use mpi
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: val(n)
integer, optional, dimension(:), intent(in) :: dims
logical, optional, intent(in) :: side
! local variables
integer comm,root_x(cart_ndims),root,ierr
logical l1
root_x=0
comm=get_cart_comm(dims)
call mpi_cart_rank(comm,root_x,root,ierr)
l1=.true.
if (present(side).and.present(dims)) then
  if (side.and..not.cart_side(dims)) l1=.false.
endif
if (l1) call mpi_bcast(val,n,MPI_DOUBLE_PRECISION,root,comm,ierr)
return
end subroutine


  
  


end module