module mod_mpi_grid

! if true then print debug information
logical debug
data debug/.false./
! number of processors
integer nproc
data nproc/1/
! index of current processor
integer iproc
data iproc/0/
! number of processors in a grid
integer mpi_grid_nproc
! number of grid dimensions
integer mpi_grid_nd
! size of each grid dimension
integer, allocatable :: mpi_grid_size(:)
! position of the current process in the grid
integer, allocatable :: mpi_grid_x(:)
! total number of grid communicators (2^number of dimensions)
!integer mpi_grid_nc
! grid communicators
integer, allocatable :: mpi_grid_comm(:)
! mpi reduce operations
integer op_sum
integer op_min
integer op_max


interface mpi_grid_bcast
  module procedure mpi_grid_bcast_d,mpi_grid_bcast_z, &
    mpi_grid_bcast_i,mpi_grid_bcast_l
end interface

interface mpi_grid_reduce
  module procedure mpi_grid_reduce_i,mpi_grid_reduce_d, &
    mpi_grid_reduce_z,mpi_grid_reduce_f,mpi_grid_reduce_i2,&
    mpi_grid_reduce_c,mpi_grid_reduce_i8
end interface

interface mpi_grid_hash
  module procedure mpi_grid_hash_z,mpi_grid_hash_d
end interface

interface mpi_grid_send
  module procedure mpi_grid_send_z,mpi_grid_send_i
end interface

interface mpi_grid_recieve
  module procedure mpi_grid_recieve_z,mpi_grid_recieve_i
end interface

contains

!--------------------------------!
!      mpi_world_initialize      !
!--------------------------------!
subroutine mpi_world_initialize
#ifdef _MPI_
use mpi
#endif
implicit none
#ifdef _MPI_
integer ierr
call mpi_init(ierr)
call mpi_comm_size(MPI_COMM_WORLD,nproc,ierr)
call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierr)
op_sum=MPI_SUM
op_min=MPI_MIN
op_max=MPI_MAX
call mpi_grid_initialize((/nproc/))
#endif
return
end subroutine

!------------------------------!
!      mpi_world_finalize      !
!------------------------------!
subroutine mpi_world_finalize
#ifdef _MPI_
use mpi
implicit none
integer ierr
call mpi_barrier(MPI_COMM_WORLD,ierr)
call mpi_finalize(ierr)
#endif
return
end subroutine

!-----------------------------!
!      mpi_world_barrier      !
!-----------------------------!
subroutine mpi_world_barrier
#ifdef _MPI_
use mpi
implicit none
integer ierr
call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif
return
end subroutine

!-------------------------------!
!      mpi_grid_initialize      !
!-------------------------------!
subroutine mpi_grid_initialize(mpi_grid_size_,debug_)
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
integer, dimension(:), intent(in) :: mpi_grid_size_
logical, optional, intent(in) :: debug_
! local variables
logical, allocatable :: l1(:)
integer ierr,i,j,nc,i1

if (present(debug_)) debug=debug_

call mpi_grid_finalize

! get number of dimensions
mpi_grid_nd=size(mpi_grid_size_)
if (allocated(mpi_grid_size)) deallocate(mpi_grid_size)
allocate(mpi_grid_size(mpi_grid_nd))
if (allocated(mpi_grid_x)) deallocate(mpi_grid_x)
allocate(mpi_grid_x(mpi_grid_nd))

#ifdef _MPI_
mpi_grid_size=mpi_grid_size_
! get number of processors in the grid
mpi_grid_nproc=1
do i=1,mpi_grid_nd
  mpi_grid_nproc=mpi_grid_nproc*mpi_grid_size(i)
enddo
if (mpi_grid_nproc.gt.nproc) then
  write(*,*)
  write(*,'("Error(mpi_grid_initialize): not enough processors to build a grid")')
  write(*,'("  number of processors : ",I8)')nproc
  write(*,'("  grid size : ",10I8)')mpi_grid_size
  call pstop
endif
! get number of communicators
nc=2**mpi_grid_nd-1
if (allocated(mpi_grid_comm)) deallocate(mpi_grid_comm)
allocate(mpi_grid_comm(0:nc))
mpi_grid_comm=MPI_COMM_NULL
mpi_grid_x=-1

allocate(l1(mpi_grid_nd))
l1=.false.

! create mpi grid
call mpi_cart_create(MPI_COMM_WORLD,mpi_grid_nd,mpi_grid_size,l1, &
  .false.,mpi_grid_comm(0),ierr)

if (mpi_grid_in()) then
! get the coordinates of the current processor
  call mpi_cart_get(mpi_grid_comm(0),mpi_grid_nd,mpi_grid_size,l1, &
    mpi_grid_x,ierr)
  if (debug.and.mpi_grid_root()) then
    write(*,*)
    write(*,'("[mpi_grid_initialize] number of grid dimensions : ",I2)')mpi_grid_nd
    write(*,'("[mpi_grid_initialize] dimension sizes : ",10I8)')mpi_grid_size   
    write(*,'("[mpi_grid_initialize] number of communicators : ",I2)')nc
  endif
! get all possible communicators
!   for example, for 3D grid we have 7 possibilities:
!     001,010,100,011,101,110,111
!   we don't have 000 (null communicator along no direction)
!   communicator 111 and global group communicator mpi_grid_comm(0)   
!   behave identically
  do i=1,nc
    l1=.false.
    i1=i
    do j=1,mpi_grid_nd
      if (mod(i1,2).eq.1) l1(j)=.true.
      i1=i1/2
    enddo
    call mpi_cart_sub(mpi_grid_comm(0),l1,mpi_grid_comm(i),ierr)
    if (debug.and.mpi_grid_root()) then
      write(*,'("[mpi_grid_initialize] index of communicator : ",I3,&
        &", communicator directions : ",10L2)')i,l1
    endif
  enddo
  mpi_grid_nproc=1
  do i=1,mpi_grid_nd
    mpi_grid_nproc=mpi_grid_nproc*mpi_grid_size(i)
  enddo
endif !mpi_grid_in
#else
! get number of dimensions
mpi_grid_size=1
mpi_grid_x=0
if (debug) then
  write(*,*)
  write(*,'("[mpi_grid_initialize] serial mode")')
endif
#endif
return
end subroutine

!-----------------------------!
!      mpi_grid_finalize      !
!-----------------------------!
subroutine mpi_grid_finalize
#ifdef _MPI_
use mpi
#endif
implicit none
#ifdef _MPI_
integer i,ierr
if (mpi_grid_in()) then
  do i=0,2**mpi_grid_nd-1
    call mpi_comm_free(mpi_grid_comm(i),ierr)
  enddo
endif
if (allocated(mpi_grid_comm)) deallocate(mpi_grid_comm)
#endif
if (allocated(mpi_grid_size)) deallocate(mpi_grid_size)
if (allocated(mpi_grid_x)) deallocate(mpi_grid_x)
return
end subroutine

!-----------------------------!
!      mpi_grid_get_comm      !
!-----------------------------!
integer function mpi_grid_get_comm(dims)
implicit none
! arguments
integer, optional, dimension(:), intent(in) :: dims
! local variables
integer i,j
#ifdef _MPI_
j=0
if (present(dims)) then
  do i=1,size(dims)
    j=j+2**(dims(i)-1)
  enddo
endif
if (debug.and.mpi_grid_root(dims).and.mpi_grid_side(dims)) then
  write(*,*)
  if (present(dims)) then
    write(*,'("[mpi_grid_get_comm] communication dimensions : ",10I3)')dims
  endif
  write(*,'("[mpi_grid_get_comm] index of communicator : ",I3)')j
endif
mpi_grid_get_comm=mpi_grid_comm(j)
#else
mpi_grid_get_comm=-1
#endif
return
end function


!----------------------------!
!      mpi_grid_barrier      !
!----------------------------!
subroutine mpi_grid_barrier(dims,side)
implicit none
integer, optional, dimension(:), intent(in) :: dims
logical, optional, intent(in) :: side
integer comm,ierr
logical side_,l1
#ifdef _MPI_
side_=.false.
if (present(side)) side_=side
l1=.true.
if (side_.and..not.mpi_grid_side(dims)) l1=.false.
comm=mpi_grid_get_comm(dims)
if (l1) call mpi_barrier(comm,ierr)
#endif
return
end subroutine

!-----------------------!
!      mpi_grid_in      !
!-----------------------!
logical function mpi_grid_in()
#ifdef _MPI_
use mpi
#endif
implicit none
#ifdef _MPI_
if (allocated(mpi_grid_comm)) then
  mpi_grid_in=(mpi_grid_comm(0).ne.MPI_COMM_NULL)
else
  mpi_grid_in=.false.
endif
#else
mpi_grid_in=.true.
#endif
return
end function

!-------------------------!
!      mpi_grid_root      !
!-------------------------!
logical function mpi_grid_root(dims)
implicit none
! arguments
integer, optional, dimension(:), intent(in) :: dims
! local variables
integer i
logical l1
#ifdef _MPI_
if (.not.present(dims)) then
  if (all(mpi_grid_x.eq.0)) then
    mpi_grid_root=.true.
  else
    mpi_grid_root=.false.
  endif
else
  l1=.true.
  do i=1,size(dims)
    if (mpi_grid_x(dims(i)).ne.0) l1=.false.
  enddo
  mpi_grid_root=l1
endif
#else
mpi_grid_root=.true.
#endif
return
end function

!-------------------------!
!      mpi_grid_side      !
!-------------------------!
! Example for 2D grid 3x4 (zero is top-left corner):
! side processors for dim=1 (dimension of size 3) are processors with 
!  second zero coordinate
!  x o o o
!  x o o o
!  x o o o  
! side processors for dim=2 (dimension of size 4) are processors with 
!  first zero coordinate
!  x x x x  
!  o o o o
!  o o o o  
logical function mpi_grid_side(dims)
implicit none
! arguments
integer, optional, dimension(:), intent(in) :: dims
! local variables
integer dims1(mpi_grid_nd)
integer i
logical l1
dims1=1
if (present(dims)) then
  do i=1,size(dims)
    dims1(dims(i))=0
  enddo
endif
l1=.true.
do i=1,mpi_grid_nd
  if (dims1(i).eq.1.and.mpi_grid_x(i).ne.0) l1=.false.
enddo
mpi_grid_side=l1
return
end function

function ortdims(dims)
implicit none
integer, intent(in) :: dims(:)
integer ortdims(mpi_grid_nd-size(dims))
integer i,j
integer f(mpi_grid_nd)
if (size(dims).eq.mpi_grid_nd) return
f=1
do j=1,size(dims)
  f(dims(j))=0
enddo
i=0
do j=1,mpi_grid_nd
  if (f(j).eq.1) then
    i=i+1
    ortdims(i)=j
  endif
enddo
return
end function

subroutine mpi_grid_bcast_common(n_,dims_,side_,root_,lbcast_,&
  length_,comm_,rootid_)
#ifdef _MPI_
use mpi
#endif
implicit none
integer, optional, intent(in) :: n_
integer, optional, dimension(:), intent(in) :: dims_
logical, optional, intent(in) :: side_
integer, optional, dimension(:), intent(in) :: root_
logical, intent(out) :: lbcast_
integer, intent(out) :: length_
integer, intent(out) :: comm_
integer, intent(out) :: rootid_
integer root_x(mpi_grid_nd),ierr,i
! check if a broadcast is necessary
if (present(dims_)) then
  lbcast_=.false.
  do i=1,size(dims_)
    if (mpi_grid_size(dims_(i)).ne.1) lbcast_=.true.
  enddo
else
  lbcast_=.true.
endif
if (.not.lbcast_) return
! check if only side processors does a broadcast
if (lbcast_) then
  if (present(side_).and.present(dims_)) then
    if (side_.and..not.mpi_grid_side(dims_)) lbcast_=.false.
  endif
endif
! length of array
length_=1
if (present(n_)) length_=n_
! get communicator
comm_=mpi_grid_get_comm(dims_)
! get root id
if (present(root_)) then
  root_x=-1
  root_x(1:size(root_))=root_
else
  root_x=0
endif
#ifdef _MPI_
call mpi_cart_rank(comm_,root_x,rootid_,ierr)
#else
rootid_=0
#endif
end subroutine


!----------------------------!
!      mpi_grid_bcast_d      !
!----------------------------!
subroutine mpi_grid_bcast_d(val,n,dims,side,root)
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
real(8), intent(in) :: val
integer, optional, intent(in) :: n
integer, optional, dimension(:), intent(in) :: dims
logical, optional, intent(in) :: side
integer, optional, dimension(:), intent(in) :: root
! local variables
integer comm,rootid,ierr,length
logical lbcast
#ifdef _MPI_
call mpi_grid_bcast_common(n_=n,dims_=dims,side_=side,root_=root,lbcast_=lbcast,&
  length_=length,comm_=comm,rootid_=rootid)
if (.not.lbcast) return
call mpi_bcast(val,length,MPI_DOUBLE_PRECISION,rootid,comm,ierr)
#endif
return
end subroutine

!----------------------------!
!      mpi_grid_bcast_z      !
!----------------------------!
subroutine mpi_grid_bcast_z(val,n,dims,side,root)
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
complex(8), intent(in) :: val
integer, optional, intent(in) :: n
integer, optional, dimension(:), intent(in) :: dims
logical, optional, intent(in) :: side
integer, optional, dimension(:), intent(in) :: root
! local variables
integer comm,rootid,ierr,length
logical lbcast
#ifdef _MPI_
call mpi_grid_bcast_common(n_=n,dims_=dims,side_=side,root_=root,lbcast_=lbcast,&
  length_=length,comm_=comm,rootid_=rootid)
if (.not.lbcast) return
call mpi_bcast(val,length,MPI_DOUBLE_COMPLEX,rootid,comm,ierr)
#endif
return
end subroutine

!----------------------------!
!      mpi_grid_bcast_i      !
!----------------------------!
subroutine mpi_grid_bcast_i(val,n,dims,side,root)
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
integer, intent(in) :: val
integer, optional, intent(in) :: n
integer, optional, dimension(:), intent(in) :: dims
logical, optional, intent(in) :: side
integer, optional, dimension(:), intent(in) :: root
! local variables
integer comm,rootid,ierr,length
logical lbcast
#ifdef _MPI_
call mpi_grid_bcast_common(n_=n,dims_=dims,side_=side,root_=root,lbcast_=lbcast,&
  length_=length,comm_=comm,rootid_=rootid)
if (.not.lbcast) return
call mpi_bcast(val,length,MPI_INTEGER,rootid,comm,ierr)
#endif
return
end subroutine

!----------------------------!
!      mpi_grid_bcast_l      !
!----------------------------!
subroutine mpi_grid_bcast_l(val,n,dims,side,root)
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
logical, intent(in) :: val
integer, optional, intent(in) :: n
integer, optional, dimension(:), intent(in) :: dims
logical, optional, intent(in) :: side
integer, optional, dimension(:), intent(in) :: root
! local variables
integer comm,rootid,ierr,length
logical lbcast
#ifdef _MPI_
call mpi_grid_bcast_common(n_=n,dims_=dims,side_=side,root_=root,lbcast_=lbcast,&
  length_=length,comm_=comm,rootid_=rootid)
if (.not.lbcast) return
call mpi_bcast(val,length,MPI_LOGICAL,rootid,comm,ierr)
#endif
return
end subroutine

subroutine mpi_grid_reduce_common(n_,dims_,side_,all_,op_,root_,lreduce_,&
  lallreduce_,length_,reduceop_,comm_,rootid_)
#ifdef _MPI_
use mpi
implicit none
integer, optional, intent(in) :: n_
integer, optional, dimension(:), intent(in) :: dims_
logical, optional, intent(in) :: side_
logical, optional, intent(in) :: all_
integer, optional, intent(in) :: op_
integer, optional, dimension(:), intent(in) :: root_
logical, intent(out) :: lreduce_
logical, intent(out) :: lallreduce_
integer, intent(out) :: length_
integer, intent(out) :: reduceop_
integer, intent(out) :: comm_
integer, intent(out) :: rootid_
integer root_x(mpi_grid_nd),ierr,i
! check if a reduction is necessary
if (present(dims_)) then
  lreduce_=.false.
  do i=1,size(dims_)
    if (mpi_grid_size(dims_(i)).ne.1) lreduce_=.true.
  enddo
else
  lreduce_=.true.
endif
! check if only side processors does a reduction
if (lreduce_) then
  if (present(side_).and.present(dims_)) then
    if (side_.and..not.mpi_grid_side(dims_)) lreduce_=.false.
  endif
endif
! check for allreduce
lallreduce_=.false.
if (present(all_)) then
  if (all_) lallreduce_=.true.
endif
! length of array
length_=1
if (present(n_)) length_=n_
! reduction operation
reduceop_=op_sum
if (present(op_)) then
  reduceop_=op_
endif
! get communicator
comm_=mpi_grid_get_comm(dims_)
! get root id
if (present(root_)) then
  root_x(1:size(root_))=root_
else
  root_x=0
endif
call mpi_cart_rank(comm_,root_x,rootid_,ierr)
if (debug) then
  write(*,'("[mpi_grid_reduce_common] lreduce_ : ",L1)')lreduce_
  write(*,'("[mpi_grid_reduce_common] lallreduce_ : ",L1)')lallreduce_
  write(*,'("[mpi_grid_reduce_common] length_ : ",I6)')length_
  write(*,'("[mpi_grid_reduce_common] reduceop_ : ",I6)')reduceop_  
  write(*,'("[mpi_grid_reduce_common] comm_ : ",I6)')comm_ 
  write(*,'("[mpi_grid_reduce_common] rootid_ : ",I6)')rootid_  
endif
#endif
end subroutine


!-----------------------------!
!      mpi_grid_reduce_i      !
!-----------------------------!
subroutine mpi_grid_reduce_i(val,n,dims,side,all,op,root)
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
integer, intent(inout) :: val
integer, optional, intent(in) :: n
integer, optional, dimension(:), intent(in) :: dims
logical, optional, intent(in) :: side
logical, optional, intent(in) :: all
integer, optional, intent(in) :: op
integer, optional, dimension(:), intent(in) :: root
! local variables
integer comm,rootid,ierr,sz
logical lreduce,lallreduce
integer length,reduceop
integer, allocatable :: tmp(:)
#ifdef _MPI_
call mpi_grid_reduce_common(n_=n,dims_=dims,side_=side,all_=all,op_=op,&
  root_=root,lreduce_=lreduce,lallreduce_=lallreduce,length_=length,&
  reduceop_=reduceop,comm_=comm,rootid_=rootid)
sz=sizeof(val)
if (.not.lreduce) return
allocate(tmp(length))
if (lallreduce) then
  call mpi_allreduce(val,tmp,length,MPI_INTEGER,reduceop,comm,ierr)
else
  call mpi_reduce(val,tmp,length,MPI_INTEGER,reduceop,rootid,comm,ierr)
endif
call memcopy(tmp,val,length*sz)
deallocate(tmp)
#endif
return
end subroutine

!-----------------------------!
!      mpi_grid_reduce_i2     !
!-----------------------------!
subroutine mpi_grid_reduce_i2(val,n,dims,side,all,op,root)
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
integer(2), intent(inout) :: val
integer, optional, intent(in) :: n
integer, optional, dimension(:), intent(in) :: dims
logical, optional, intent(in) :: side
logical, optional, intent(in) :: all
integer, optional, intent(in) :: op
integer, optional, dimension(:), intent(in) :: root
! local variables
integer comm,rootid,ierr,sz
logical lreduce,lallreduce
integer length,reduceop
integer(2), allocatable :: tmp(:)
#ifdef _MPI_
call mpi_grid_reduce_common(n_=n,dims_=dims,side_=side,all_=all,op_=op,&
  root_=root,lreduce_=lreduce,lallreduce_=lallreduce,length_=length,&
  reduceop_=reduceop,comm_=comm,rootid_=rootid)
sz=sizeof(val)
if (.not.lreduce) return
allocate(tmp(length))
if (lallreduce) then
  call mpi_allreduce(val,tmp,length,MPI_INTEGER2,reduceop,comm,ierr)
else
  call mpi_reduce(val,tmp,length,MPI_INTEGER2,reduceop,rootid,comm,ierr)
endif
call memcopy(tmp,val,length*sz)
deallocate(tmp)
#endif
return
end subroutine

!-----------------------------!
!      mpi_grid_reduce_i8     !
!-----------------------------!
subroutine mpi_grid_reduce_i8(val,n,dims,side,all,op,root)
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
integer(8), intent(inout) :: val
integer, optional, intent(in) :: n
integer, optional, dimension(:), intent(in) :: dims
logical, optional, intent(in) :: side
logical, optional, intent(in) :: all
integer, optional, intent(in) :: op
integer, optional, dimension(:), intent(in) :: root
! local variables
integer comm,rootid,ierr,sz
logical lreduce,lallreduce
integer length,reduceop
integer(8), allocatable :: tmp(:)
#ifdef _MPI_
call mpi_grid_reduce_common(n_=n,dims_=dims,side_=side,all_=all,op_=op,&
  root_=root,lreduce_=lreduce,lallreduce_=lallreduce,length_=length,&
  reduceop_=reduceop,comm_=comm,rootid_=rootid)
sz=sizeof(val)
if (.not.lreduce) return
allocate(tmp(length))
if (lallreduce) then
  call mpi_allreduce(val,tmp,length,MPI_INTEGER8,reduceop,comm,ierr)
else
  call mpi_reduce(val,tmp,length,MPI_INTEGER8,reduceop,rootid,comm,ierr)
endif
call memcopy(tmp,val,length*sz)
deallocate(tmp)
#endif
return
end subroutine

!-----------------------------!
!      mpi_grid_reduce_d      !
!-----------------------------!
subroutine mpi_grid_reduce_d(val,n,dims,side,all,op,root)
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
real(8), intent(inout) :: val
integer, optional, intent(in) :: n
integer, optional, dimension(:), intent(in) :: dims
logical, optional, intent(in) :: side
logical, optional, intent(in) :: all
integer, optional, intent(in) :: op
integer, optional, dimension(:), intent(in) :: root
! local variables
integer comm,rootid,ierr,sz
logical lreduce,lallreduce
integer length,reduceop
real(8), allocatable :: tmp(:)
#ifdef _MPI_
call mpi_grid_reduce_common(n_=n,dims_=dims,side_=side,all_=all,op_=op,&
  root_=root,lreduce_=lreduce,lallreduce_=lallreduce,length_=length,&
  reduceop_=reduceop,comm_=comm,rootid_=rootid)
sz=sizeof(val)
if (.not.lreduce) return
allocate(tmp(length))
if (lallreduce) then
  call mpi_allreduce(val,tmp,length,MPI_DOUBLE_PRECISION,reduceop,comm,ierr)
else
  call mpi_reduce(val,tmp,length,MPI_DOUBLE_PRECISION,reduceop,rootid,comm,ierr)
endif
call memcopy(tmp,val,length*sz)
deallocate(tmp)
#endif
return
end subroutine

!-----------------------------!
!      mpi_grid_reduce_z      !
!-----------------------------!
subroutine mpi_grid_reduce_z(val,n,dims,side,all,op,root)
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
complex(8), intent(inout) :: val
integer, optional, intent(in) :: n
integer, optional, dimension(:), intent(in) :: dims
logical, optional, intent(in) :: side
logical, optional, intent(in) :: all
integer, optional, intent(in) :: op
integer, optional, dimension(:), intent(in) :: root
! local variables
integer comm,rootid,ierr,sz
logical lreduce,lallreduce
integer length,reduceop
complex(8), allocatable :: tmp(:)
#ifdef _MPI_
if (debug) then
  write(*,'("[mpi_grid_reduce_z] start")')
endif
call mpi_grid_reduce_common(n_=n,dims_=dims,side_=side,all_=all,op_=op,&
  root_=root,lreduce_=lreduce,lallreduce_=lallreduce,length_=length,&
  reduceop_=reduceop,comm_=comm,rootid_=rootid)
sz=sizeof(val)
if (.not.lreduce) return
allocate(tmp(length))
if (lallreduce) then
  call mpi_allreduce(val,tmp,length,MPI_DOUBLE_COMPLEX,reduceop,comm,ierr)
else
  call mpi_reduce(val,tmp,length,MPI_DOUBLE_COMPLEX,reduceop,rootid,comm,ierr)
endif
call memcopy(tmp,val,length*sz)
deallocate(tmp)
#endif
return
end subroutine

!-----------------------------!
!      mpi_grid_reduce_c      !
!-----------------------------!
subroutine mpi_grid_reduce_c(val,n,dims,side,all,op,root)
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
complex(4), intent(inout) :: val
integer, optional, intent(in) :: n
integer, optional, dimension(:), intent(in) :: dims
logical, optional, intent(in) :: side
logical, optional, intent(in) :: all
integer, optional, intent(in) :: op
integer, optional, dimension(:), intent(in) :: root
! local variables
integer comm,rootid,ierr,sz
logical lreduce,lallreduce
integer length,reduceop
complex(4), allocatable :: tmp(:)
#ifdef _MPI_
if (debug) then
  write(*,'("[mpi_grid_reduce_c] start")')
endif
call mpi_grid_reduce_common(n_=n,dims_=dims,side_=side,all_=all,op_=op,&
  root_=root,lreduce_=lreduce,lallreduce_=lallreduce,length_=length,&
  reduceop_=reduceop,comm_=comm,rootid_=rootid)
sz=sizeof(val)
if (.not.lreduce) return
allocate(tmp(length))
if (lallreduce) then
  call mpi_allreduce(val,tmp,length,MPI_COMPLEX,reduceop,comm,ierr)
else
  call mpi_reduce(val,tmp,length,MPI_COMPLEX,reduceop,rootid,comm,ierr)
endif
call memcopy(tmp,val,length*sz)
deallocate(tmp)
#endif
return
end subroutine

!-----------------------------!
!      mpi_grid_reduce_f      !
!-----------------------------!
subroutine mpi_grid_reduce_f(val,n,dims,side,all,op,root)
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
real(4), intent(inout) :: val
integer, optional, intent(in) :: n
integer, optional, dimension(:), intent(in) :: dims
logical, optional, intent(in) :: side
logical, optional, intent(in) :: all
integer, optional, intent(in) :: op
integer, optional, dimension(:), intent(in) :: root
! local variables
integer comm,rootid,ierr,sz
logical lreduce,lallreduce
integer length,reduceop
real(4), allocatable :: tmp(:)
#ifdef _MPI_
call mpi_grid_reduce_common(n_=n,dims_=dims,side_=side,all_=all,op_=op,&
  root_=root,lreduce_=lreduce,lallreduce_=lallreduce,length_=length,&
  reduceop_=reduceop,comm_=comm,rootid_=rootid)
sz=sizeof(val)
if (.not.lreduce) return
allocate(tmp(length))
if (lallreduce) then
  call mpi_allreduce(val,tmp,length,MPI_REAL,reduceop,comm,ierr)
else
  call mpi_reduce(val,tmp,length,MPI_REAL,reduceop,rootid,comm,ierr)
endif
call memcopy(tmp,val,length*sz)
deallocate(tmp)
#endif
return
end subroutine

!---------------------------!
!      mpi_grid_hash_z      !
!---------------------------!
subroutine mpi_grid_hash_z(val,n,d,ierr)
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
complex(8), intent(in) :: val
integer, intent(in) :: n
integer, intent(in) :: d
integer, intent(out) :: ierr
! local variables
integer, allocatable :: tmp(:)
integer, external :: hash
integer i
#ifdef _MPI_
ierr=0
allocate(tmp(0:mpi_grid_size(d)-1))
tmp=0
tmp(mpi_grid_x(d))=hash(val,16*n)
call mpi_grid_reduce(tmp(0),mpi_grid_size(d),dims=(/d/),all=.true.)
if (mpi_grid_x(d).eq.0) then
  do i=0,mpi_grid_size(d)-1
    if (tmp(i).ne.tmp(0)) ierr=1
  enddo
else
  if (tmp(mpi_grid_x(d)).ne.tmp(0)) ierr=1
endif
deallocate(tmp)
#endif
return
end subroutine

!---------------------------!
!      mpi_grid_hash_d      !
!---------------------------!
subroutine mpi_grid_hash_d(val,n,d,ierr)
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
real(8), intent(in) :: val
integer, intent(in) :: n
integer, intent(in) :: d
integer, intent(out) :: ierr
! local variables
integer, allocatable :: tmp(:)
integer, external :: hash
integer i
#ifdef _MPI_
ierr=0
allocate(tmp(0:mpi_grid_size(d)-1))
tmp=0
tmp(mpi_grid_x(d))=hash(val,8*n)
call mpi_grid_reduce(tmp(0),mpi_grid_size(d),dims=(/d/),all=.true.)
if (mpi_grid_x(d).eq.0) then
  do i=0,mpi_grid_size(d)-1
    if (tmp(i).ne.tmp(0)) ierr=1
  enddo
else
  if (tmp(mpi_grid_x(d)).ne.tmp(0)) ierr=1
endif
deallocate(tmp)
#endif
return
end subroutine


subroutine mpi_grid_msg(sname,msg)
implicit none
character*(*), intent(in) :: sname
character*(*), intent(in) :: msg
character*100 sx
character*256 sout
character*20 c1
integer i
sx="x"
do i=1,mpi_grid_nd
  write(c1,'(I10)')mpi_grid_x(i)
  sx=trim(adjustl(sx))//":"//trim(adjustl(c1))
enddo
sout="["//trim(adjustl(sname))//" "//trim(adjustl(sx))//"] "//trim(adjustl(msg))
write(*,'(A)')trim(adjustl(sout))
return
end subroutine


!------------------------!
!      mpi_grid_map      !
!------------------------!
! Examples of use
!  get size of local block for current processor
!    nkptloc=mpi_grid_map(nkp,dim_k)
!  get global index by local index of current processor
!    ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
!  get local index of k-point for current processor
!    ikloc=mpi_grid_map(nkp,dim_k,glob=ik)
!  get local index of k-point for processor with coordinate h
!    ikloc=mpi_grid_map(nkp,dim_k,x=h,glob=ik)
integer function mpi_grid_map(length,idim,x,loc,glob,offs)
implicit none
! arguments
integer, intent(in) :: length
integer, intent(in) :: idim
integer, optional, intent(inout) :: x
integer, optional, intent(in) :: loc
integer, optional, intent(in) :: glob
integer, optional, intent(out) :: offs
! local variables
integer idx0_,size_,x_
integer i

i=0
if (present(loc)) i=i+1
if (present(glob)) i=i+1
if (present(offs)) i=i+1

if (i.gt.1) then
  write(*,'("Error(cart_map): more than one optional argument is presented")')
  call pstop
endif

if (present(glob)) then
  call idxloc(length,mpi_grid_size(idim),glob,x_,idx0_)
  if (present(x)) x=x_-1
  mpi_grid_map=idx0_
  return
endif

if (present(x)) then
  x_=x
else
  x_=mpi_grid_x(idim)
endif

if (present(loc)) then
  call idxglob(length,mpi_grid_size(idim),x_+1,loc,idx0_)
  mpi_grid_map=idx0_
  return
endif

call idxofs(length,mpi_grid_size(idim),x_+1,idx0_,size_)
if (present(offs)) offs=idx0_
mpi_grid_map=size_
return
end function

!---------------------------!
!      mpi_grid_send_z      !
!---------------------------!
subroutine mpi_grid_send_z(val,n,dims,dest,tag)
#ifdef _MPI_
use mpi
#endif
implicit none
complex(8), intent(in) :: val
integer, intent(in) :: n
integer, dimension(:), intent(in) :: dims
integer, dimension(:), intent(in) :: dest
integer, intent(in) :: tag
! local variables
integer comm,dest_rank,req,ierr
if (debug) then
  write(*,*)'[mpi_grid_send_z] mpi_grid_x:',mpi_grid_x
  write(*,*)'[mpi_grid_send_z] n=',n
  write(*,*)'[mpi_grid_send_z] dims=',dims
  write(*,*)'[mpi_grid_send_z] dest=',dest
  write(*,*)'[mpi_grid_send_z] tag=',tag
endif
#ifdef _MPI_
comm=mpi_grid_get_comm(dims)
call mpi_cart_rank(comm,dest,dest_rank,ierr) 
call mpi_isend(val,n,MPI_DOUBLE_COMPLEX,dest_rank,tag,comm,req,ierr)
#endif
return
end subroutine

!---------------------------!
!      mpi_grid_send_i      !
!---------------------------!
subroutine mpi_grid_send_i(val,n,dims,dest,tag)
#ifdef _MPI_
use mpi
#endif
implicit none
integer, intent(in) :: val
integer, intent(in) :: n
integer, dimension(:), intent(in) :: dims
integer, dimension(:), intent(in) :: dest
integer, intent(in) :: tag
! local variables
integer comm,dest_rank,req,ierr
#ifdef _MPI_
comm=mpi_grid_get_comm(dims)
call mpi_cart_rank(comm,dest,dest_rank,ierr) 
call mpi_isend(val,n,MPI_INTEGER,dest_rank,tag,comm,req,ierr)
#endif
return
end subroutine

!------------------------------!
!      mpi_grid_recieve_z      !
!------------------------------!
subroutine mpi_grid_recieve_z(val,n,dims,src,tag)
#ifdef _MPI_
use mpi
#endif
implicit none
complex(8), intent(out) :: val
integer, intent(in) :: n
integer, dimension(:), intent(in) :: dims
integer, dimension(:), intent(in) :: src
integer, intent(in) :: tag
! local variables
#ifdef _MPI_
integer comm,src_rank,ierr
integer stat(MPI_STATUS_SIZE)
#endif
if (debug) then
  write(*,*)'[mpi_grid_recieve_z] mpi_grid_x:',mpi_grid_x
  write(*,*)'[mpi_grid_recieve_z] n=',n
  write(*,*)'[mpi_grid_recieve_z] dims=',dims
  write(*,*)'[mpi_grid_recieve_z] src=',src
  write(*,*)'[mpi_grid_recieve_z] tag=',tag
endif
#ifdef _MPI_
comm=mpi_grid_get_comm(dims)
call mpi_cart_rank(comm,src,src_rank,ierr)
call mpi_recv(val,n,MPI_DOUBLE_COMPLEX,src_rank,tag,comm,stat,ierr)
#endif
return
end subroutine

!------------------------------!
!      mpi_grid_recieve_i      !
!------------------------------!
subroutine mpi_grid_recieve_i(val,n,dims,src,tag)
#ifdef _MPI_
use mpi
#endif
implicit none
integer, intent(out) :: val
integer, intent(in) :: n
integer, dimension(:), intent(in) :: dims
integer, dimension(:), intent(in) :: src
integer, intent(in) :: tag
! local variables
#ifdef _MPI_
integer comm,src_rank,ierr
integer stat(MPI_STATUS_SIZE)
comm=mpi_grid_get_comm(dims)
call mpi_cart_rank(comm,src,src_rank,ierr)
call mpi_recv(val,n,MPI_INTEGER,src_rank,tag,comm,stat,ierr)
#endif
return
end subroutine





! partition something of size "length" to blocks 
!   (in) nblocks is number of blocks
!   (in) iblock is the index of block (from 1 to nblocks)
!   (out) idx0 is the offset in global block
!   (out) blocksize is the size of local block
! if length.eq.0 -> nothing to partition
! if iblock outside the boundaries -> return with error 
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
if (length.lt.0) lerr=.true.
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
if (length.eq.0) then
  idx0=-1
  blocksize=0
  return
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

!-----------------!
!      pstop      !
!-----------------!
subroutine pstop(ierr_)
#ifdef _MPI_
use mpi
#endif
implicit none
integer, optional, intent(in) :: ierr_
integer ierr

if (present(ierr_)) then
  ierr=ierr_
else
  ierr=-1
endif

write(*,'("STOP execution")')
write(*,'("  error code : ",I8)')ierr
write(*,'("  global index of processor : ",I8)')iproc
if (allocated(mpi_grid_x)) &
  write(*,'("  coordinates of processor : ",10I8)')mpi_grid_x
#ifdef _MPI_
call mpi_abort(MPI_COMM_WORLD,ierr,ierr)
call mpi_finalize(ierr)
#endif
stop
return
end subroutine

!-----------------!
!      bstop      !
!-----------------!
subroutine bstop(ierr_)
#ifdef _MPI_
use mpi
#endif
implicit none
integer, optional, intent(in) :: ierr_
integer ierr

if (present(ierr_)) then
  ierr=ierr_
else
  ierr=-1
endif

write(*,'("STOP execution")')
write(*,'("  error code : ",I8)')ierr
write(*,'("  global index of processor : ",I8)')iproc
if (allocated(mpi_grid_x)) &
  write(*,'("  coordinates of processor : ",10I8)')mpi_grid_x
#ifdef _MPI_
call mpi_grid_barrier()
call mpi_finalize(ierr)
#endif
stop
return
end subroutine

end module

subroutine memcopy(src,dest,size)
implicit none
character, intent(in) :: src(*)
character, intent(out) :: dest(*)
integer, intent(in) :: size
dest(1:size)=src(1:size)
return
end




