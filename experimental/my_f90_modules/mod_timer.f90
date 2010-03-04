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
!integer, allocatable :: cart_map_size(:)

interface cart_bcast
  module procedure d_cart_bcast,z_cart_bcast,l_cart_bcast
end interface

interface cart_reduce
  module procedure i_cart_reduce,d_cart_reduce
end interface

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

subroutine cart_barrier(dims)
implicit none
integer, optional, dimension(:), intent(in) :: dims
integer comm,ierr
comm=get_cart_comm(dims)
call mpi_barrier(comm,ierr)
return
end subroutine

subroutine pstop(ierr_)
use mpi
implicit none
integer, optional, intent(in) :: ierr_
integer ierr
write(*,'("STOP parallel execution")')
write(*,'("  global index of processor : ",I6)')iproc
if (allocated(cart_x)) &
  write(*,'("  Cartesian coordinates of processor : ",10I4)')cart_x
if (present(ierr_)) then
  ierr=ierr_
else
  ierr=-1
endif
call mpi_abort(MPI_COMM_WORLD,ierr,ierr)
call mpi_finalize(ierr)
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
real(8), intent(in) :: val
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

subroutine l_cart_bcast(val,n,dims,side)
use mpi
implicit none
! arguments
logical, intent(in) :: val
integer, optional, intent(in) :: n
integer, optional, dimension(:), intent(in) :: dims
logical, optional, intent(in) :: side
! local variables
integer comm,root_x(cart_ndims),root,ierr,n_
logical l1
n_=1
if (present(n)) n_=n
root_x=0
comm=get_cart_comm(dims)
call mpi_cart_rank(comm,root_x,root,ierr)
l1=.true.
if (present(side).and.present(dims)) then
  if (side.and..not.cart_side(dims)) l1=.false.
endif
if (l1) call mpi_bcast(val,n_,MPI_LOGICAL,root,comm,ierr)
return
end subroutine

!--------------------------------!
!      complex(8) broadcast      !
!--------------------------------!
subroutine z_cart_bcast(val,n,dims,side)
use mpi
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: val
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
if (l1) call mpi_bcast(val,n,MPI_DOUBLE_COMPLEX,root,comm,ierr)
return
end subroutine

!-----------------------------!
!      integer reduction      !
!-----------------------------!
subroutine i_cart_reduce(val,n,dims,side,all,op)
use mpi
implicit none
! arguments
integer, intent(in) :: n
integer, intent(inout) :: val
integer, optional, dimension(:), intent(in) :: dims
logical, optional, intent(in) :: side
logical, optional, intent(in) :: all
integer, optional, intent(in) :: op
! local variables
integer comm,root_x(cart_ndims),root,ierr,sz
logical all_,l1
integer op_
integer, allocatable :: tmp(:)

sz=sizeof(val)

all_=.false.
if (present(all)) then
  if (all) all_=.true.
endif
op_=MPI_SUM
if (present(op)) then
  op_=op
endif

comm=get_cart_comm(dims)
root_x=0
call mpi_cart_rank(comm,root_x,root,ierr)
l1=.true.
if (present(side).and.present(dims)) then
  if (side.and..not.cart_side(dims)) l1=.false.
endif
if (l1) then
  if (all_) then
    allocate(tmp(n))
    call mpi_allreduce(val,tmp,n,MPI_INTEGER,op_,comm,ierr)
    call memcopy(tmp,val,n*sz)
    deallocate(tmp)
  else
    if (cart_root(dims)) allocate(tmp(n))
    call mpi_reduce(val,tmp,n,MPI_INTEGER,op_,root,comm,ierr)
    if (cart_root(dims)) then
      call memcopy(tmp,val,n*sz)
      deallocate(tmp)
    endif
  endif
endif
return
end subroutine

!-----------------------------!
!      real(8) reduction      !
!-----------------------------!
subroutine d_cart_reduce(val,n,dims,side,all,op)
use mpi
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(inout) :: val
integer, optional, dimension(:), intent(in) :: dims
logical, optional, intent(in) :: side
logical, optional, intent(in) :: all
integer, optional, intent(in) :: op
! local variables
integer comm,root_x(cart_ndims),root,ierr,sz
logical all_,l1
integer op_
real(8), allocatable :: tmp(:)

sz=sizeof(val)

all_=.false.
if (present(all)) then
  if (all) all_=.true.
endif
op_=MPI_SUM
if (present(op)) then
  op_=op
endif

comm=get_cart_comm(dims)
root_x=0
call mpi_cart_rank(comm,root_x,root,ierr)
l1=.true.
if (present(side).and.present(dims)) then
  if (side.and..not.cart_side(dims)) l1=.false.
endif
if (l1) then
  if (all_) then
    allocate(tmp(n))
    call mpi_allreduce(val,tmp,n,MPI_DOUBLE_PRECISION,op_,comm,ierr)
    call memcopy(tmp,val,n*sz)
    deallocate(tmp)
  else
    if (cart_root(dims)) allocate(tmp(n))
    call mpi_reduce(val,tmp,n,MPI_DOUBLE_PRECISION,op_,root,comm,ierr)
    if (cart_root(dims)) then
      call memcopy(tmp,val,n*sz)
      deallocate(tmp)
    endif
  endif
endif
return
end subroutine


integer function cart_map(length,idim,x,loc,glob,offs)
implicit none
! arguments
integer, intent(in) :: length
integer, intent(in) :: idim
integer, optional, intent(inout) :: x
integer, optional, intent(inout) :: loc
integer, optional, intent(inout) :: glob
integer, optional, intent(out) :: offs
! local variables
integer idx0_,size_,x_
integer i

!if (.not.cart_ingrid()) then
!  cart_map=0
!  return
!endif

i=0
if (present(loc)) i=i+1
if (present(glob)) i=i+1
if (present(offs)) i=i+1

if (i.gt.1) then
  write(*,'("Error(cart_map) : more than one optional argument is presented")')
  call pstop
endif

if (present(glob)) then
  call idxloc(length,cart_dims(idim),glob,x_,idx0_)
  if (present(x)) x=x_
  cart_map=idx0_
  return
endif

if (present(x)) then
  x_=x
else
  x_=cart_x(idim)
endif

if (present(loc)) then
  call idxglob(length,cart_dims(idim),x_+1,loc,idx0_)
  cart_map=idx0_
  return
endif

call idxofs(length,cart_dims(idim),x_+1,idx0_,size_)
if (present(offs)) offs=idx0_
cart_map=size_
return
end function

subroutine idxofs(length,nblocks,iblock,idx0,blocksize)
implicit none
integer, intent(in) :: length
integer, intent(in) :: nblocks
integer, intent(in) :: iblock
integer, intent(out) :: idx0
integer, intent(out) :: blocksize

integer n1,n2
logical lerr

lerr=.false.
if (length.lt.1) lerr=.true.
if (nblocks.lt.1) lerr=.true.
if (iblock.lt.1.or.iblock.gt.nblocks) lerr=.true.
if (lerr) then
  write(*,*)
  write(*,'("Error(idxofs) : wrong input arguments")')
  write(*,'("  length : ",I8)')length
  write(*,'("  nblocks : ",I8)')nblocks
  write(*,'("  iblock : ",I8)')iblock
  write(*,*)
  call pstop
endif
n1=length/nblocks
n2=mod(length,nblocks)
if (n1.eq.0) then
  if (iblock.le.n2) then
    idx0=iblock-1
    blocksize=1
  else
    idx0=-1
    blocksize=0
  endif
else
  if (iblock.le.n2) then
    idx0=(n1+1)*(iblock-1)
    blocksize=n1+1
  else
    idx0=(n1+1)*n2+n1*(iblock-n2-1)
    blocksize=n1
  endif
endif
return
end subroutine

subroutine idxloc(length,nblocks,idxg,iblock,idxl)
implicit none
integer, intent(in) :: length
integer, intent(in) :: nblocks
integer, intent(in) :: idxg
integer, intent(out) :: iblock
integer, intent(out) :: idxl

integer n1,n2
logical lerr

lerr=.false.
if (length.lt.1) lerr=.true.
if (nblocks.lt.1) lerr=.true.
if (idxg.lt.1.or.idxg.gt.length) lerr=.true.
if (lerr) then
  write(*,*)
  write(*,'("Error(idxloc) : wrong input arguments")')
  write(*,'("  length : ",I8)')length
  write(*,'("  nblocks : ",I8)')nblocks
  write(*,'("  idxg : ",I8)')idxg
  write(*,*)
  call pstop
endif
n1=length/nblocks
n2=mod(length,nblocks)
if (idxg.le.(n1+1)*n2) then
  iblock=(idxg-1)/(n1+1)+1
  idxl=mod(idxg-1,n1+1)+1
else
  iblock=n2+(idxg-(n1+1)*n2-1)/n1+1
  idxl=mod(idxg-(n1+1)*n2-1,n1)+1
endif
return
end subroutine

subroutine idxglob(length,nblocks,iblock,idxl,idxg)
implicit none
integer, intent(in) :: length
integer, intent(in) :: nblocks
integer, intent(in) :: iblock
integer, intent(in) :: idxl
integer, intent(out) :: idxg

integer n1,n2
logical lerr

lerr=.false.
if (length.lt.1) lerr=.true.
if (nblocks.lt.1) lerr=.true.
if (iblock.lt.1.or.iblock.gt.nblocks) lerr=.true.
if (idxl.lt.1) lerr=.true.
if (lerr) then
  write(*,*)
  write(*,'("Error(idxglob) : wrong input arguments")')
  write(*,'("  length : ",I8)')length
  write(*,'("  nblocks : ",I8)')nblocks
  write(*,'("  iblock : ",I8)')iblock
  write(*,'("  idxl : ",I8)')idxl
  write(*,*)
  call pstop
endif
n1=length/nblocks
n2=mod(length,nblocks)
lerr=.false.
if (iblock.le.n2.and.(idxl.gt.(n1+1))) lerr=.true.
if (iblock.gt.n2.and.(idxl.gt.n1)) lerr=.true.
if (lerr) then
  write(*,*)
  write(*,'("Error(idxloc) : local index out of boundary")')
  write(*,'("  iblock : ",I8)')iblock
  write(*,'("  idxl : ",I8)')idxl
  write(*,*)
  call pstop
endif  
if (iblock.le.n2) then
  idxg=(iblock-1)*(n1+1)+idxl
else
  idxg=n2*(n1+1)+(iblock-n2-1)*n1+idxl
endif
return
end subroutine

  
  


end module