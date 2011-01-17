!> @brief Fortran module to work with MPI Cartesian grids
!! @details MPI grid is useful when you have several independent
!! variables which can be parallelized. MPI specification provides
!! basic API to work with the grids. This module masks MPI calls and
!! adds more flexibility to the basic MPI routines. 
!! @author Anton Kozhevnikov
!! @date 2009
module mod_mpi_grid

!> print debug information if set to .true.
logical, private :: mpi_grid_debug
data mpi_grid_debug/.false./
!> number of processes allocated for the job
integer nproc
data nproc/1/
!> index of the current process (starting from 0)
integer iproc
data iproc/0/
!> number of processes in the grid
integer mpi_grid_nproc
!> maximum number of grid dimensions
integer, parameter :: ndmax=5
!> number of grid dimensions
integer mpi_grid_nd
!> number of grid communicators
integer mpi_grid_nc
!> size of each grid dimension
integer, private, allocatable :: mpi_grid_size(:)
!> coordinates [0;mpi_grid_size(i)-1) of the current process in the grid
integer, private, allocatable :: mpi_grid_x(:)
!> grid communicators
integer, private, allocatable :: mpi_grid_comm(:)
!> mpi reduction operations: sum
integer op_sum
!> mpi reduction operations: min
integer op_min
!> mpi reduction operations: max
integer op_max

@python ftypes=["integer(2)","integer(4)","integer(8)","real(4)","real(8)","complex(4)","complex(8)","logical"]
@python fsuffixes=["_i2","_i4","_i8","_f","_d","_c","_z","_l"]
@python fmpitypes=["MPI_INTEGER2","MPI_INTEGER","MPI_INTEGER8","MPI_REAL","MPI_DOUBLE_PRECISION","MPI_COMPLEX","MPI_DOUBLE_COMPLEX","MPI_LOGICAL"]
@python ntypes=8

!> @brief Interface to MPI reduction operation
!> @param inpb input buffer
!> @param n (optional) buffer size
!> @param outb (optional) output buffer
!> @param dims (optional) communication dimensions 
!> @param side (optional) if .true. then only side processes do a reduction
!> @param all (optional) if .true. then mpi_allreduce is executed
!> @param op (optional) reduction operation
!> @param root (optional) coordinates of the "root" process, which gathers 
!! the full reduced buffer
interface mpi_grid_reduce
@template begin
@template variable fsuffix
@python for i in range(ntypes-1): fsuffix=fsuffixes[i];
  module procedure mpi_grid_reduce#fsuffix
@template end
end interface

interface mpi_grid_bcast
@template begin
@template variable fsuffix
@python for i in range(ntypes): fsuffix=fsuffixes[i];
  module procedure mpi_grid_bcast#fsuffix
@template end
end interface

interface mpi_grid_hash
@template begin
@template variable fsuffix
@python for i in range(ntypes-1): fsuffix=fsuffixes[i];
  module procedure mpi_grid_hash#fsuffix
@template end
end interface

interface mpi_grid_send
@template begin
@template variable fsuffix
@python for i in range(ntypes): fsuffix=fsuffixes[i];
  module procedure mpi_grid_send#fsuffix
@template end
end interface

interface mpi_grid_recieve
@template begin
@template variable fsuffix
@python for i in range(ntypes): fsuffix=fsuffixes[i];
  module procedure mpi_grid_recieve#fsuffix
@template end
end interface

! public API
public mpi_initialize
public mpi_world_initialize
public mpi_world_finalize
public mpi_world_barrier
public mpi_grid_initialize
public mpi_grid_finalize
public mpi_grid_in
public mpi_grid_root
public mpi_grid_barrier
public mpi_grid_map
public bstop
public ortdims
public mpi_grid_dim_size
public mpi_grid_dim_pos
public mpi_grid_reduce
public mpi_grid_bcast
public mpi_grid_hash
public mpi_grid_send
public mpi_grid_recieve
! private functions, used internally by the module
private convert_dims_to_internal
private mpi_grid_root_internal
private mpi_grid_side_internal
private mpi_grid_get_comm_internal
private mpi_grid_rank_internal
private mpi_grid_reduce_common
@template begin
@template variable fsuffix
@python for i in range(ntypes-1): fsuffix=fsuffixes[i];
private mpi_grid_reduce#fsuffix
@template end
private mpi_grid_bcast_common
@template begin
@template variable fsuffix
@python for i in range(ntypes): fsuffix=fsuffixes[i];
private mpi_grid_bcast#fsuffix
@template end
@template begin
@template variable fsuffix
@python for i in range(ntypes-1): fsuffix=fsuffixes[i];
private mpi_grid_hash#fsuffix
@template end
@template begin
@template variable fsuffix
@python for i in range(ntypes): fsuffix=fsuffixes[i];
private mpi_grid_send#fsuffix
@template end
@template begin
@template variable fsuffix
@python for i in range(ntypes): fsuffix=fsuffixes[i];
private mpi_grid_recieve#fsuffix
@template end
private partition_index
private local_index
private global_index

contains

!> @brief Initialize MPI library.
!> @details This subroutine calls mpi_init
subroutine mpi_initialize
#ifdef _MPI_
use mpi
#endif
implicit none
#ifdef _MPI_
integer ierr
call mpi_init(ierr)
#endif
return
end subroutine

!> @brief Initialize MPI world.
!> @details Get the total number of MPI processes, index of the current 
!! process and MPI reduction operations
subroutine mpi_world_initialize
#ifdef _MPI_
use mpi
#endif
implicit none
#ifdef _MPI_
integer ierr
call mpi_comm_size(MPI_COMM_WORLD,nproc,ierr)
call mpi_comm_rank(MPI_COMM_WORLD,iproc,ierr)
op_sum=MPI_SUM
op_min=MPI_MIN
op_max=MPI_MAX
#endif
return
end subroutine

!> @brief Finalize MPI library usage.
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

!> @brief Global barrier.
subroutine mpi_world_barrier
#ifdef _MPI_
use mpi
implicit none
integer ierr
call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif
return
end subroutine

!> @brief Initialize Cartesian MPI grid.
!> @param grid_size size of the array determines the dimensionality of the
!! MPI grid and elements - length of each dimension
!> @param grid_debug (optional) .true. if debug information should be printed
subroutine mpi_grid_initialize(grid_size,grid_debug)
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
integer, dimension(:), intent(in) :: grid_size
logical, optional, intent(in) :: grid_debug
! local variables
logical, allocatable :: l1(:)
integer ierr,i,j,nc,i1

if (present(grid_debug)) mpi_grid_debug=grid_debug

call mpi_grid_finalize

! get number of dimensions
mpi_grid_nd=size(grid_size)
if (mpi_grid_nd.gt.ndmax) then
  write(*,'("Error(mpi_grid_initialize): too many grid dimensions")')
  write(*,'("  input grid size : ",10I8)')grid_size
  call pstop
endif
allocate(mpi_grid_size(mpi_grid_nd))
allocate(mpi_grid_x(mpi_grid_nd))

#ifdef _MPI_
mpi_grid_size=grid_size
! get number of processes in the grid
mpi_grid_nproc=1
do i=1,mpi_grid_nd
  mpi_grid_nproc=mpi_grid_nproc*mpi_grid_size(i)
enddo
if (mpi_grid_nproc.gt.nproc) then
  write(*,*)
  write(*,'("Error(mpi_grid_initialize): not enough processes to build a grid")')
  write(*,'("  number of processes : ",I8)')nproc
  write(*,'("  grid size : ",10I8)')mpi_grid_size
  call pstop
endif
! get number of communicators
nc=2**mpi_grid_nd
allocate(mpi_grid_comm(nc))
mpi_grid_comm=MPI_COMM_NULL
mpi_grid_x=-1
mpi_grid_nc=nc

allocate(l1(mpi_grid_nd))
l1=.false.
! create mpi grid
call mpi_cart_create(MPI_COMM_WORLD,mpi_grid_nd,mpi_grid_size,l1, &
  .false.,mpi_grid_comm(nc),ierr)

if (mpi_grid_in()) then
! get the coordinates of the current process
  call mpi_cart_get(mpi_grid_comm(nc),mpi_grid_nd,mpi_grid_size,l1, &
    mpi_grid_x,ierr)
  if (mpi_grid_debug.and.mpi_grid_root()) then
    write(*,*)
    write(*,'("[mpi_grid_initialize] number of grid dimensions : ",I2)')mpi_grid_nd
    write(*,'("[mpi_grid_initialize] dimension sizes : ",10I8)')mpi_grid_size   
    write(*,'("[mpi_grid_initialize] number of communicators : ",I2)')nc
  endif
! get all possible communicators
!   for example, for 3D grid we have 7 possibilities:
!     001,010,100,011,101,110,111
!   we don't have 000 (null communicator along no direction)
!   communicator 111 and global group communicator mpi_grid_comm(nc)   
!   behave identically
  do i=1,nc-1
    l1=.false.
    i1=i
    do j=1,mpi_grid_nd
      if (mod(i1,2).eq.1) l1(j)=.true.
      i1=i1/2
    enddo
    call mpi_cart_sub(mpi_grid_comm(nc),l1,mpi_grid_comm(i),ierr)
    if (mpi_grid_debug.and.mpi_grid_root()) then
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
if (mpi_grid_debug) then
  write(*,*)
  write(*,'("[mpi_grid_initialize] serial mode")')
endif
#endif
return
end subroutine

!> @brief Finalize Cartesian MPI grid
subroutine mpi_grid_finalize
#ifdef _MPI_
use mpi
#endif
implicit none
#ifdef _MPI_
integer i,ierr
if (mpi_grid_in()) then
  do i=1,mpi_grid_nc
    call mpi_comm_free(mpi_grid_comm(i),ierr)
  enddo
endif
if (allocated(mpi_grid_comm)) deallocate(mpi_grid_comm)
#endif
if (allocated(mpi_grid_size)) deallocate(mpi_grid_size)
if (allocated(mpi_grid_x)) deallocate(mpi_grid_x)
return
end subroutine

!> @brief .true. if process is in the grid
logical function mpi_grid_in()
#ifdef _MPI_
use mpi
#endif
implicit none
#ifdef _MPI_
if (allocated(mpi_grid_comm)) then
  mpi_grid_in=(mpi_grid_comm(mpi_grid_nc).ne.MPI_COMM_NULL)
else
  mpi_grid_in=.false.
endif
#else
mpi_grid_in=.true.
#endif
return
end function

!> @brief .true. if this process is the root for a given list of dimensions
!> @param dims (optional) variable-length array with dimensions; if not provided,
!! the whole grid is assumed
!> @details Suppose we have a 2-dimensional 3x4 grid (3 rows, 4 columns) : 
!> @verbatim
!> o o o o
!> o o o o
!> o o o o
!> @endverbatim
!> Example usage:
!> @code
!> if (mpi_grid_root((/1,3/))) then
!>   write(*,*)'root processes of (1,3) sub-grid have coordinates : ',mpi_grid_x
!> end if
!> @endcode
logical function mpi_grid_root(dims)
implicit none
! arguments
integer, optional, dimension(:), intent(in) :: dims
! local variables
integer idims(0:ndmax)
idims=convert_dims_to_internal(dims)
mpi_grid_root=mpi_grid_root_internal(idims)
return
end function

!> @brief Barrier for a given sub-grid
!> @param dims (optional) variable-length array with dimensions; if not provided,
!! the whole grid is assumed
!> @param side (optional) if .true. then only side processes do the barrier
subroutine mpi_grid_barrier(dims,side)
implicit none
integer, optional, dimension(:), intent(in) :: dims
logical, optional, intent(in) :: side
integer comm,ierr
logical side_,l1
integer idims(0:ndmax)
#ifdef _MPI_
idims=convert_dims_to_internal(dims)
side_=.false.
if (present(side)) side_=side
l1=.true.
if (side_.and..not.mpi_grid_side_internal(idims)) l1=.false.
if (idims(0).ne.0) then
  comm=mpi_grid_get_comm_internal(idims)
  if (l1) call mpi_barrier(comm,ierr)
endif
#endif
return
end subroutine

!> @brief .true. if processes are at the side of the sub-grid
!> @param dims (optional) variable-length array with dimensions; if not provided,
!! the whole grid is assumed
!> @details Suppose we have a 2-dimensional 3x4 grid (3 rows, 4 columns) : 
!> @verbatim
!> o o o o
!> o o o o
!> o o o o
!> @endverbatim
!> side processes for dim=1 (dimension of size 3) are the processes with 
!! second zero coordinate
!> @verbatim
!> x o o o
!> x o o o
!> x o o o  
!> @endverbatim
!> side processes for dim=2 (dimension of size 4) are processes with 
!! first zero coordinate
!> @verbatim
!> x x x x  
!> o o o o
!> o o o o  
!> @endverbatim
!> Example usage:
!> @code
!> if (mpi_grid_side((/1,3/))) then
!>   write(*,*)'side processes of (1,3) sub-grid have coordinates : ',mpi_grid_x
!> end if
!> @endcode
logical function mpi_grid_side(dims)
implicit none
integer, optional, dimension(:), intent(in) :: dims
integer idims(0:ndmax)
idims=convert_dims_to_internal(dims)
mpi_grid_side=mpi_grid_side_internal(idims)
return
end function

!> @brief return orthogonal subset of dimensions
!> @param dims (optional) variable-length array with dimensions; if not provided,
!! the whole grid is assumed
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


subroutine mpi_grid_msg(sname,msg)
implicit none
character*(*), intent(in) :: sname
character*(*), intent(in) :: msg
character*100 sx
character*1000 sout
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



! private subroutines and functions


!> @brief Private function of the module
!> @details Covert variable-length array with dimensions to a fixed-size array
!> internally used in the module
!> @param dims variable-length array with dimensions
function convert_dims_to_internal(dims)
implicit none
integer, optional, dimension(:), intent(in) :: dims
integer, dimension(0:ndmax) :: convert_dims_to_internal
integer n,i
if (present(dims)) then
  convert_dims_to_internal=0
  n=0
  do i=1,size(dims)
    if (dims(i).le.mpi_grid_nd) then
      n=n+1
      if (n.gt.ndmax) then
        write(*,'("Error(convert_dims_to_internal): too many grid dimensions")')
        write(*,'("  dims : ",10I8)')dims
        call pstop
      endif
      convert_dims_to_internal(0)=n
      convert_dims_to_internal(n)=dims(i)
    endif
  enddo
else
  convert_dims_to_internal(0)=mpi_grid_nd
  do i=1,mpi_grid_nd
    convert_dims_to_internal(i)=i
  enddo
endif
return
end function

!> @brief Private function of the module
!> @details .true. if this process is the root for a given list of dimensions
!> @param idims fixed-size array with dimensions
logical function mpi_grid_root_internal(idims)
implicit none
! arguments
integer, intent(in) :: idims(0:ndmax)
! local variables
integer i
logical l1
#ifdef _MPI_
l1=.true.
do i=1,idims(0)
  if (mpi_grid_x(idims(i)).ne.0) l1=.false.
enddo
mpi_grid_root_internal=l1
#else
mpi_grid_root_internal=.true.
#endif
return
end function

!> @brief Private function of the module
!> @details return communicator for a given sub-grid
!> @param idims fixed-size array with dimensions
integer function mpi_grid_get_comm_internal(idims)
implicit none
! arguments
integer, intent(in) :: idims(0:ndmax)
! local variables
integer i,j
#ifdef _MPI_
if (idims(0).eq.0) then
  write(*,'("Error(mpi_grid_get_comm_internal): no dimensions")')
  call pstop
endif
j=0
do i=1,idims(0)
  j=j+2**(idims(i)-1)
enddo
if (mpi_grid_debug.and.mpi_grid_root_internal(idims).and.&
    mpi_grid_side_internal(idims)) then
  write(*,*)
  write(*,'("[mpi_grid_get_comm_internal] communication dimensions : ",10I3)')&
    idims(1:idims(0))
  write(*,'("[mpi_grid_get_comm_internal] index of communicator : ",I3)')j
endif
mpi_grid_get_comm_internal=mpi_grid_comm(j)
#else
mpi_grid_get_comm_internal=-1
#endif
return
end function

!> @brief Private function of the module
!> @details .true. if processes are at the side of the sub-grid
!> @param idims fixed-size array with dimensions
logical function mpi_grid_side_internal(idims)
implicit none
! arguments
integer, intent(in) :: idims(0:ndmax)
! local variables
integer dims1(mpi_grid_nd)
integer i
logical l1
dims1=1
do i=1,idims(0)
  dims1(idims(i))=0
enddo
l1=.true.
do i=1,mpi_grid_nd
  if (dims1(i).eq.1.and.mpi_grid_x(i).ne.0) l1=.false.
enddo
mpi_grid_side_internal=l1
return
end function

integer function mpi_grid_rank_internal(comm,ndims,x)
#ifdef _MPI_
use mpi
implicit none
integer, intent(in) :: comm
integer, intent(in) :: ndims
integer, optional, dimension(:), intent(in) :: x
integer x1(mpi_grid_nd),k,rank,ierr,i
! default position
x1=0
if (present(x)) then
  k=min(ndims,size(x))
  do i=1,k
    x1(i)=x(i)
  enddo
endif
do i=1,mpi_grid_nd
  if (x1(i).ge.mpi_grid_size(i)) then
    write(*,'("Error(mpi_grid_rank_internal): process coordinates are &
      &out of range")')
    write(*,'("  coordinates : ",10I8)')x1
    write(*,'("  mpi_grid_size : ",10I8)')mpi_grid_size
    call pstop
  endif
enddo
call mpi_cart_rank(comm,x1,rank,ierr)
mpi_grid_rank_internal=rank
#endif
return
end function

subroutine mpi_grid_bcast_common(n_,dims_,side_,root_,lbcast_,&
  length_,comm_,rootid_)
#ifdef _MPI_
use mpi
implicit none
integer, optional, intent(in) :: n_
integer, optional, dimension(:), intent(in) :: dims_
logical, optional, intent(in) :: side_
integer, optional, dimension(:), intent(in) :: root_
logical, intent(out) :: lbcast_
integer, intent(out) :: length_
integer, intent(out) :: comm_
integer, intent(out) :: rootid_
integer i,idims(0:ndmax)
!
idims=convert_dims_to_internal(dims_)
! check if broadcast is necessary
lbcast_=.false.
do i=1,idims(0)
  if (mpi_grid_size(idims(i)).ne.1) lbcast_=.true.
enddo
! check if only side processes do a broadcast
if (lbcast_) then
  if (present(side_)) then
    if (side_.and..not.mpi_grid_side_internal(idims)) lbcast_=.false.
  endif
endif
if (.not.lbcast_) return
! length of array
length_=1
if (present(n_)) length_=n_
! get communicator
comm_=mpi_grid_get_comm_internal(idims)
! get root id
rootid_=mpi_grid_rank_internal(comm_,idims(0),root_)
#endif
end subroutine

@template begin
@template variable fsuffix
@template variable ftype
@template variable fmpitype
@python for i in range(ntypes): fsuffix=fsuffixes[i]; ftype=ftypes[i]; fmpitype=fmpitypes[i];
subroutine mpi_grid_bcast#fsuffix(val,n,dims,side,root)
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
#ftype, intent(in) :: val
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
call mpi_bcast(val,length,#fmpitype,rootid,comm,ierr)
#endif
return
end subroutine

@template end

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
integer idims(0:ndmax)
!
idims=convert_dims_to_internal(dims_)
! check if a reduction is necessary
!lreduce_=.false.
!do i=1,idims(0)
!  if (mpi_grid_size(idims(i)).ne.1) lreduce_=.true.
!enddo
lreduce_=.true.
! check if only side processes do a reduction
if (lreduce_) then
  if (present(side_)) then
    if (side_.and..not.mpi_grid_side_internal(idims)) lreduce_=.false.
  endif
endif
if (.not.lreduce_) return
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
comm_=mpi_grid_get_comm_internal(idims)
! get root id
rootid_=mpi_grid_rank_internal(comm_,idims(0),root_)
if (mpi_grid_debug) then
  write(*,'("[mpi_grid_reduce_common] lreduce_ : ",L1)')lreduce_
  write(*,'("[mpi_grid_reduce_common] lallreduce_ : ",L1)')lallreduce_
  write(*,'("[mpi_grid_reduce_common] length_ : ",I6)')length_
  write(*,'("[mpi_grid_reduce_common] reduceop_ : ",I6)')reduceop_  
  write(*,'("[mpi_grid_reduce_common] comm_ : ",I6)')comm_ 
  write(*,'("[mpi_grid_reduce_common] rootid_ : ",I6)')rootid_  
endif
#endif
end subroutine

@template begin
@template variable fsuffix
@template variable ftype
@template variable fmpitype
@python for i in range(ntypes-1): fsuffix=fsuffixes[i]; ftype=ftypes[i]; fmpitype=fmpitypes[i];
!> @brief Private function of the module
!> @details #ftype implementation of mpi_grid_reduce
!> @param idims fixed-size array with dimensions
subroutine mpi_grid_reduce#fsuffix(inpb,n,outb,dims,side,all,op,root)
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
#ftype, intent(inout) :: inpb
integer, optional, intent(in) :: n
#ftype, optional, intent(out) :: outb
integer, optional, dimension(:), intent(in) :: dims
logical, optional, intent(in) :: side
logical, optional, intent(in) :: all
integer, optional, intent(in) :: op
integer, optional, dimension(:), intent(in) :: root
! local variables
integer comm,rootid,ierr,sz
logical lreduce,lallreduce
integer length,reduceop
#ftype, allocatable :: tmp(:)
#ifdef _MPI_
call mpi_grid_reduce_common(n_=n,dims_=dims,side_=side,all_=all,op_=op,&
  root_=root,lreduce_=lreduce,lallreduce_=lallreduce,length_=length,&
  reduceop_=reduceop,comm_=comm,rootid_=rootid)
if (.not.lreduce) return
! do a reduction
if (present(outb)) then
  if (lallreduce) then
    call mpi_allreduce(inpb,outb,length,#fmpitype,reduceop,comm,ierr)
  else
    call mpi_reduce(inpb,outb,length,#fmpitype,reduceop,rootid,comm,ierr)
  endif
else
  allocate(tmp(length))
  if (lallreduce) then
    call mpi_allreduce(inpb,tmp,length,#fmpitype,reduceop,comm,ierr)
  else
    call mpi_reduce(inpb,tmp,length,#fmpitype,reduceop,rootid,comm,ierr)
  endif
  sz=sizeof(inpb)
  call memcopy(tmp,inpb,length*sz)
  deallocate(tmp)
endif
#endif
return
end subroutine

@template end

@template begin
@template variable fsuffix
@template variable ftype
@template variable fmpitype
@python for i in range(ntypes-1): fsuffix=fsuffixes[i]; ftype=ftypes[i]; fmpitype=fmpitypes[i];
subroutine mpi_grid_hash#fsuffix(val,n,d,ierr)
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
#ftype, intent(in) :: val
integer, intent(in) :: n
integer, intent(in) :: d
integer, intent(out) :: ierr
! local variables
integer, allocatable :: tmp(:)
integer, external :: hash
integer i,sz
#ifdef _MPI_
ierr=0
allocate(tmp(0:mpi_grid_size(d)-1))
tmp=0
sz=sizeof(val)
tmp(mpi_grid_x(d))=hash(val,n*sz)
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

@template end

@template begin
@template variable fsuffix
@template variable ftype
@template variable fmpitype
@python for i in range(ntypes): fsuffix=fsuffixes[i]; ftype=ftypes[i]; fmpitype=fmpitypes[i];
subroutine mpi_grid_send#fsuffix(val,n,dims,dest,tag)
#ifdef _MPI_
use mpi
#endif
implicit none
#ftype, intent(in) :: val
integer, intent(in) :: n
integer, dimension(:), intent(in) :: dims
integer, dimension(:), intent(in) :: dest
integer, intent(in) :: tag
! local variables
integer comm,dest_rank,req,ierr
integer idims(0:ndmax)
idims=convert_dims_to_internal(dims)
if (mpi_grid_debug) then
  write(*,*)'[mpi_grid_send#fsuffix] mpi_grid_x:',mpi_grid_x
  write(*,*)'[mpi_grid_send#fsuffix] n=',n
  write(*,*)'[mpi_grid_send#fsuffix] dims=',dims
  write(*,*)'[mpi_grid_send#fsuffix] dest=',dest
  write(*,*)'[mpi_grid_send#fsuffix] tag=',tag
endif
#ifdef _MPI_
comm=mpi_grid_get_comm_internal(idims)
dest_rank=mpi_grid_rank_internal(comm,idims(0),dest)
call mpi_isend(val,n,#fmpitype,dest_rank,tag,comm,req,ierr)
#endif
return
end subroutine

@template end

@template begin
@template variable fsuffix
@template variable ftype
@template variable fmpitype
@python for i in range(ntypes): fsuffix=fsuffixes[i]; ftype=ftypes[i]; fmpitype=fmpitypes[i];
subroutine mpi_grid_recieve#fsuffix(val,n,dims,src,tag)
#ifdef _MPI_
use mpi
#endif
implicit none
#ftype, intent(out) :: val
integer, intent(in) :: n
integer, dimension(:), intent(in) :: dims
integer, dimension(:), intent(in) :: src
integer, intent(in) :: tag
! local variables
#ifdef _MPI_
integer comm,src_rank,ierr
integer stat(MPI_STATUS_SIZE)
integer idims(0:ndmax)
idims=convert_dims_to_internal(dims)
#endif
if (mpi_grid_debug) then
  write(*,*)'[mpi_grid_recieve#fsuffix] mpi_grid_x:',mpi_grid_x
  write(*,*)'[mpi_grid_recieve#fsuffix] n=',n
  write(*,*)'[mpi_grid_recieve#fsuffix] dims=',dims
  write(*,*)'[mpi_grid_recieve#fsuffix] src=',src
  write(*,*)'[mpi_grid_recieve#fsuffix] tag=',tag
endif
#ifdef _MPI_
comm=mpi_grid_get_comm_internal(idims)
src_rank=mpi_grid_rank_internal(comm,idims(0),src)
call mpi_recv(val,n,#fmpitype,src_rank,tag,comm,stat,ierr)
#endif
return
end subroutine

@template end


integer function mpi_grid_dim_size(idim)
implicit none
! arguments
integer, intent(in) :: idim
!
if (idim.gt.mpi_grid_nd) then
  mpi_grid_dim_size=1
else
  mpi_grid_dim_size=mpi_grid_size(idim)
endif
return
end function

integer function mpi_grid_dim_pos(idim)
implicit none
! arguments
integer, intent(in) :: idim
!
if (idim.gt.mpi_grid_nd) then
  mpi_grid_dim_pos=0
else
  mpi_grid_dim_pos=mpi_grid_x(idim)
endif
return
end function

!------------------------!
!      mpi_grid_map      !
!------------------------!
! Examples of use
!  get size of local block for current process
!    nkptloc=mpi_grid_map(nkp,dim_k)
!  get global index by local index of current process
!    ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
!  get local index of k-point for current process
!    ikloc=mpi_grid_map(nkp,dim_k,glob=ik)
!  get local index of k-point for process with coordinate h
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
integer idx0_,size_,x_,dimsz
integer i

i=0
if (present(loc)) i=i+1
if (present(glob)) i=i+1
if (present(offs)) i=i+1

if (i.gt.1) then
  write(*,'("Error(mpi_grid_map): more than one optional argument is presented")')
  call pstop
endif

dimsz=mpi_grid_dim_size(idim)

if (present(glob)) then
  call local_index(length,dimsz,glob,x_,idx0_)
  if (present(x)) x=x_-1
  mpi_grid_map=idx0_
  return
endif

if (present(x)) then
  x_=x
else
  x_=mpi_grid_dim_pos(idim)
endif

if (x_.lt.0.or.x_.ge.dimsz) then
  write(*,'("Error(mpi_grid_map): coordinate out of range")')
  write(*,'("  x_ : ",I8)')x_
  call pstop
endif

if (present(loc)) then
  call global_index(length,dimsz,x_+1,loc,idx0_)
  mpi_grid_map=idx0_
  return
endif

call partition_index(length,dimsz,x_+1,idx0_,size_)
if (present(offs)) offs=idx0_
mpi_grid_map=size_
return
end function


integer function mpi_grid_map2(length,dims,x,loc,glob,offs)
implicit none
! arguments
integer, intent(in) :: length
integer, dimension(:), intent(in) :: dims
integer, optional, dimension(:), intent(inout) :: x
integer, optional, intent(in) :: loc
integer, optional, intent(in) :: glob
integer, optional, intent(out) :: offs
! local variables
integer idx0_,size_
integer i,size1,k,j
integer, allocatable :: x_(:)

i=0
if (present(loc)) i=i+1
if (present(glob)) i=i+1
if (present(offs)) i=i+1

if (i.gt.1) then
  write(*,'("Error(mpi_grid_map2): more than one optional argument is presented")')
  call pstop
endif

! total number of processes
size1=1
do i=1,size(dims)
  size1=size1*mpi_grid_size(dims(i))
enddo

!if (present(glob)) then
!  call idxloc(length,mpi_grid_size(idim),glob,x_,idx0_)
!  if (present(x)) x=x_-1
!  mpi_grid_map=idx0_
!  return
!endif

allocate(x_(size(dims)))
if (present(x)) then
  if (size(x).ne.size(dims)) then
    write(*,'("Error(mpi_grid_map2): size(x).ne.size(dims)")')
    call pstop
  endif
  do i=1,size(dims)
    x_(i)=x(i)
  enddo
else
  do i=1,size(dims)
    x_(i)=mpi_grid_x(dims(i))
  enddo
endif 

! make a liner index from coordinates
k=1
j=1
do i=1,size(dims)
  j=j+x_(i)*k
  k=k*mpi_grid_size(dims(i))
enddo

!if (present(loc)) then
!  call idxglob(length,mpi_grid_size(idim),x_+1,loc,idx0_)
!  mpi_grid_map=idx0_
!  return
!endif

call partition_index(length,size1,j,idx0_,size_)
if (present(offs)) offs=idx0_
mpi_grid_map2=size_

deallocate(x_)
return
end function

!> @brief Partition global index
!> @param sz (in) global size of the index
!> @param np (in) number of partitions
!> @param ip (in) i-th partition
!> @param ofs (out) offset in global index
!> @param psz (out) size of i-th partition
!> @details global index (1 2 3 4 5) is partitioned into two local blocks
!! as (1 2 3) (4 5); into three block as (1 2) (3 4) (5)
!> @verbatim
!> index           | 1  2  3  4  5 |  
!> partitions      | 1   | 2   | 3 |
!> partition sizes | 2   | 2   | 1 |
!> offsets         | 0   | 2   | 4 |
!> @endverbatim
subroutine partition_index(sz,np,ip,ofs,psz)
implicit none
! arguments
integer, intent(in) :: sz
integer, intent(in) :: np
integer, intent(in) :: ip
integer, intent(out) :: ofs
integer, intent(out) :: psz
! local variables
integer n1,n2
logical lerr
! 
lerr=.false.
if (sz.lt.0) lerr=.true.
if (np.lt.1) lerr=.true.
if (ip.lt.1.or.ip.gt.np) lerr=.true.
if (lerr) then
  write(*,*)
  write(*,'("Error(partition_index) : wrong input arguments")')
  write(*,'("  sz : ",I8)')sz
  write(*,'("  np : ",I8)')np
  write(*,'("  ip : ",I8)')ip
  write(*,*)
  call pstop
endif
if (sz.eq.0) then
  ofs=-1
  psz=0
  return
endif
! minimum partition size  
n1=sz/np
! first n2 partitions get extra element
n2=mod(sz,np)
if (ip.le.n2) then
  ofs=(n1+1)*(ip-1)
  psz=n1+1
else
  ofs=(n1+1)*n2+n1*(ip-n2-1)
  psz=n1
endif
return
end subroutine

subroutine local_index(sz,np,idxg,ip,idxl)
implicit none
! arguments
integer, intent(in) :: sz
integer, intent(in) :: np
integer, intent(in) :: idxg
integer, intent(out) :: ip
integer, intent(out) :: idxl
! local variables
integer n1,n2
logical lerr
!
lerr=.false.
if (sz.lt.1) lerr=.true.
if (np.lt.1) lerr=.true.
if (idxg.lt.1.or.idxg.gt.sz) lerr=.true.
if (lerr) then
  write(*,*)
  write(*,'("Error(local_index) : wrong input arguments")')
  write(*,'("  sz : ",I8)')sz
  write(*,'("  np : ",I8)')np
  write(*,'("  idxg : ",I8)')idxg
  write(*,*)
  call pstop
endif
n1=sz/np
n2=mod(sz,np)
if (idxg.le.(n1+1)*n2) then
  ip=(idxg-1)/(n1+1)+1
  idxl=mod(idxg-1,n1+1)+1
else
  ip=n2+(idxg-(n1+1)*n2-1)/n1+1
  idxl=mod(idxg-(n1+1)*n2-1,n1)+1
endif
return
end subroutine

subroutine global_index(sz,np,ip,idxl,idxg)
implicit none
! arguments
integer, intent(in) :: sz
integer, intent(in) :: np
integer, intent(in) :: ip
integer, intent(in) :: idxl
integer, intent(out) :: idxg
! local variables
integer n1,n2
logical lerr
!
lerr=.false.
if (sz.lt.1) lerr=.true.
if (np.lt.1) lerr=.true.
if (ip.lt.1.or.ip.gt.np) lerr=.true.
if (idxl.lt.1) lerr=.true.
if (lerr) then
  write(*,*)
  write(*,'("Error(global_index) : wrong input arguments")')
  write(*,'("  sz : ",I8)')sz
  write(*,'("  np : ",I8)')np
  write(*,'("  ip : ",I8)')ip
  write(*,'("  idxl : ",I8)')idxl
  write(*,*)
  call pstop
endif
n1=sz/np
n2=mod(sz,np)
lerr=.false.
if (ip.le.n2.and.(idxl.gt.(n1+1))) lerr=.true.
if (ip.gt.n2.and.(idxl.gt.n1)) lerr=.true.
if (lerr) then
  write(*,*)
  write(*,'("Error(global_index) : local index out of boundary")')
  write(*,'("  ip : ",I8)')ip
  write(*,'("  idxl : ",I8)')idxl
  write(*,*)
  call pstop
endif  
if (ip.le.n2) then
  idxg=(ip-1)*(n1+1)+idxl
else
  idxg=n2*(n1+1)+(ip-n2-1)*n1+idxl
endif
return
end subroutine

!-----------------!
!      bstop      !
!-----------------!
subroutine bstop
#ifdef _MPI_
use mpi
#endif
implicit none
integer ierr
write(*,'("STOP execution")')
write(*,'("  global index of process : ",I8)')iproc
if (allocated(mpi_grid_x)) &
  write(*,'("  coordinates of process : ",10I8)')mpi_grid_x
#ifdef _MPI_
call mpi_grid_barrier()
call mpi_finalize(ierr)
#endif
stop
return
end subroutine

end module

!-----------------!
!      pstop      !
!-----------------!
subroutine pstop
#ifdef _MPI_
use mpi
use mod_mpi_grid
implicit none
integer ierr
write(*,'("STOP execution")')
write(*,'("  global index of process : ",I8)')iproc
call mpi_abort(MPI_COMM_WORLD,-1,ierr)
call mpi_finalize(ierr)
#endif
stop
return
end subroutine

subroutine memcopy(src,dest,size)
implicit none
character, intent(in) :: src(*)
character, intent(out) :: dest(*)
integer, intent(in) :: size
dest(1:size)=src(1:size)
return
end




