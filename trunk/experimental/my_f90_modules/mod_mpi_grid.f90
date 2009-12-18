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


!interface mpi_grid_bcast
!  module procedure d_cart_bcast,z_cart_bcast,l_cart_bcast
!end interface

!interface cart_reduce
!  module procedure i_cart_reduce,d_cart_reduce
!end interface

contains

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
#endif
return
end subroutine

subroutine mpi_world_finalize
#ifdef _MPI_
use mpi
#endif
implicit none
#ifdef _MPI_
integer ierr
call mpi_barrier(MPI_COMM_WORLD,ierr)
call mpi_finalize(ierr)
#endif
return
end subroutine

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

#ifdef _MPI_
! get number of dimensions
mpi_grid_nd=size(mpi_grid_size_)
if (allocated(mpi_grid_size)) deallocate(mpi_grid_size)
allocate(mpi_grid_size(mpi_grid_nd))
mpi_grid_size=mpi_grid_size_
! get number of communicators
nc=2**mpi_grid_nd-1
if (allocated(mpi_grid_comm)) deallocate(mpi_grid_comm)
allocate(mpi_grid_comm(0:nc))
mpi_grid_comm=MPI_COMM_NULL
if (allocated(mpi_grid_x)) deallocate(mpi_grid_x)
allocate(mpi_grid_x(mpi_grid_nd))
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
    write(*,'("[mod_mpi_grid] number of grid dimensions : ",I2)')mpi_grid_nd
    write(*,'("[mod_mpi_grid] dimension sizes : ",10I8)')mpi_grid_size   
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
      write(*,'("[mod_mpi_grid] index of communicator : ",I3,&
        &", communicator directions : ",10L2)')i,l1
    endif
  enddo
endif !mpi_grid_in
#else
! get number of dimensions
mpi_grid_nd=1
if (allocated(mpi_grid_size)) deallocate(mpi_grid_size)
allocate(mpi_grid_size(mpi_grid_nd))
mpi_grid_size=1
if (allocated(mpi_grid_x)) deallocate(mpi_grid_x)
allocate(mpi_grid_x(mpi_grid_nd))
mpi_grid_x=0
if (debug) then
  write(*,*)
  write(*,'("[mod_mpi_grid] in serial mode")')
endif
#endif
return
end subroutine

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
deallocate(mpi_grid_comm)
#endif
deallocate(mpi_grid_size)
deallocate(mpi_grid_x)
return
end subroutine

!subroutine cart_barrier(dims)
!implicit none
!integer, optional, dimension(:), intent(in) :: dims
!integer comm,ierr
!comm=get_cart_comm(dims)
!call mpi_barrier(comm,ierr)
!return
!end subroutine
!

logical function mpi_grid_in()
#ifdef _MPI_
use mpi
#endif
implicit none
#ifdef _MPI_
mpi_grid_in=(mpi_grid_comm(0).ne.MPI_COMM_NULL)
#else
mpi_grid_in=.true.
#endif
return
end function

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
!
!logical function cart_side(dims)
!implicit none
!! arguments
!integer, dimension(:), intent(in) :: dims
!! local variables
!integer dims1(cart_ndims)
!integer i
!logical l1
!dims1=1
!do i=1,size(dims)
!  dims1(dims(i))=0
!enddo
!l1=.true.
!do i=1,cart_ndims
!  if (dims1(i).eq.1.and.cart_x(i).ne.0) l1=.false.
!enddo
!cart_side=l1
!return
!end function
!
!integer function get_cart_comm(dims)
!implicit none
!! arguments
!integer, optional, dimension(:), intent(in) :: dims
!! local variables
!integer i,j
!
!if (.not.present(dims)) then
!  get_cart_comm=cart_comm(0)
!else
!  j=1
!  do i=1,size(dims)
!    j=j+2**(dims(i)-1)
!  enddo
!  get_cart_comm=cart_comm(j)
!endif
!return
!end function
!
!subroutine d_cart_bcast(val,n,dims,side)
!use mpi
!implicit none
!! arguments
!integer, intent(in) :: n
!real(8), intent(in) :: val
!integer, optional, dimension(:), intent(in) :: dims
!logical, optional, intent(in) :: side
!! local variables
!integer comm,root_x(cart_ndims),root,ierr
!logical l1
!root_x=0
!comm=get_cart_comm(dims)
!call mpi_cart_rank(comm,root_x,root,ierr)
!l1=.true.
!if (present(side).and.present(dims)) then
!  if (side.and..not.cart_side(dims)) l1=.false.
!endif
!if (l1) call mpi_bcast(val,n,MPI_DOUBLE_PRECISION,root,comm,ierr)
!return
!end subroutine
!
!subroutine l_cart_bcast(val,n,dims,side)
!use mpi
!implicit none
!! arguments
!logical, intent(in) :: val
!integer, optional, intent(in) :: n
!integer, optional, dimension(:), intent(in) :: dims
!logical, optional, intent(in) :: side
!! local variables
!integer comm,root_x(cart_ndims),root,ierr,n_
!logical l1
!n_=1
!if (present(n)) n_=n
!root_x=0
!comm=get_cart_comm(dims)
!call mpi_cart_rank(comm,root_x,root,ierr)
!l1=.true.
!if (present(side).and.present(dims)) then
!  if (side.and..not.cart_side(dims)) l1=.false.
!endif
!if (l1) call mpi_bcast(val,n_,MPI_LOGICAL,root,comm,ierr)
!return
!end subroutine
!
!!--------------------------------!
!!      complex(8) broadcast      !
!!--------------------------------!
!subroutine z_cart_bcast(val,n,dims,side)
!use mpi
!implicit none
!! arguments
!integer, intent(in) :: n
!complex(8), intent(in) :: val
!integer, optional, dimension(:), intent(in) :: dims
!logical, optional, intent(in) :: side
!! local variables
!integer comm,root_x(cart_ndims),root,ierr
!logical l1
!root_x=0
!comm=get_cart_comm(dims)
!call mpi_cart_rank(comm,root_x,root,ierr)
!l1=.true.
!if (present(side).and.present(dims)) then
!  if (side.and..not.cart_side(dims)) l1=.false.
!endif
!if (l1) call mpi_bcast(val,n,MPI_DOUBLE_COMPLEX,root,comm,ierr)
!return
!end subroutine
!
!!-----------------------------!
!!      integer reduction      !
!!-----------------------------!
!subroutine i_cart_reduce(val,n,dims,side,all,op)
!use mpi
!implicit none
!! arguments
!integer, intent(in) :: n
!integer, intent(inout) :: val
!integer, optional, dimension(:), intent(in) :: dims
!logical, optional, intent(in) :: side
!logical, optional, intent(in) :: all
!integer, optional, intent(in) :: op
!! local variables
!integer comm,root_x(cart_ndims),root,ierr,sz
!logical all_,l1
!integer op_
!integer, allocatable :: tmp(:)
!
!sz=sizeof(val)
!
!all_=.false.
!if (present(all)) then
!  if (all) all_=.true.
!endif
!op_=MPI_SUM
!if (present(op)) then
!  op_=op
!endif
!
!comm=get_cart_comm(dims)
!root_x=0
!call mpi_cart_rank(comm,root_x,root,ierr)
!l1=.true.
!if (present(side).and.present(dims)) then
!  if (side.and..not.cart_side(dims)) l1=.false.
!endif
!if (l1) then
!  if (all_) then
!    allocate(tmp(n))
!    call mpi_allreduce(val,tmp,n,MPI_INTEGER,op_,comm,ierr)
!    call memcopy(tmp,val,n*sz)
!    deallocate(tmp)
!  else
!    if (cart_root(dims)) allocate(tmp(n))
!    call mpi_reduce(val,tmp,n,MPI_INTEGER,op_,root,comm,ierr)
!    if (cart_root(dims)) then
!      call memcopy(tmp,val,n*sz)
!      deallocate(tmp)
!    endif
!  endif
!endif
!return
!end subroutine
!
!!-----------------------------!
!!      real(8) reduction      !
!!-----------------------------!
!subroutine d_cart_reduce(val,n,dims,side,all,op)
!use mpi
!implicit none
!! arguments
!integer, intent(in) :: n
!real(8), intent(inout) :: val
!integer, optional, dimension(:), intent(in) :: dims
!logical, optional, intent(in) :: side
!logical, optional, intent(in) :: all
!integer, optional, intent(in) :: op
!! local variables
!integer comm,root_x(cart_ndims),root,ierr,sz
!logical all_,l1
!integer op_
!real(8), allocatable :: tmp(:)
!
!sz=sizeof(val)
!
!all_=.false.
!if (present(all)) then
!  if (all) all_=.true.
!endif
!op_=MPI_SUM
!if (present(op)) then
!  op_=op
!endif
!
!comm=get_cart_comm(dims)
!root_x=0
!call mpi_cart_rank(comm,root_x,root,ierr)
!l1=.true.
!if (present(side).and.present(dims)) then
!  if (side.and..not.cart_side(dims)) l1=.false.
!endif
!if (l1) then
!  if (all_) then
!    allocate(tmp(n))
!    call mpi_allreduce(val,tmp,n,MPI_DOUBLE_PRECISION,op_,comm,ierr)
!    call memcopy(tmp,val,n*sz)
!    deallocate(tmp)
!  else
!    if (cart_root(dims)) allocate(tmp(n))
!    call mpi_reduce(val,tmp,n,MPI_DOUBLE_PRECISION,op_,root,comm,ierr)
!    if (cart_root(dims)) then
!      call memcopy(tmp,val,n*sz)
!      deallocate(tmp)
!    endif
!  endif
!endif
!return
!end subroutine
!
!
!integer function cart_map(length,idim,x,loc,glob,offs)
!implicit none
!! arguments
!integer, intent(in) :: length
!integer, intent(in) :: idim
!integer, optional, intent(inout) :: x
!integer, optional, intent(inout) :: loc
!integer, optional, intent(inout) :: glob
!integer, optional, intent(out) :: offs
!! local variables
!integer idx0_,size_,x_
!integer i
!
!!if (.not.cart_ingrid()) then
!!  cart_map=0
!!  return
!!endif
!
!i=0
!if (present(loc)) i=i+1
!if (present(glob)) i=i+1
!if (present(offs)) i=i+1
!
!if (i.gt.1) then
!  write(*,'("Error(cart_map) : more than one optional argument is presented")')
!  call pstop
!endif
!
!if (present(glob)) then
!  call idxloc(length,cart_dims(idim),glob,x_,idx0_)
!  if (present(x)) x=x_
!  cart_map=idx0_
!  return
!endif
!
!if (present(x)) then
!  x_=x
!else
!  x_=cart_x(idim)
!endif
!
!if (present(loc)) then
!  call idxglob(length,cart_dims(idim),x_+1,loc,idx0_)
!  cart_map=idx0_
!  return
!endif
!
!call idxofs(length,cart_dims(idim),x_+1,idx0_,size_)
!if (present(offs)) offs=idx0_
!cart_map=size_
!return
!end function
!
!subroutine idxofs(length,nblocks,iblock,idx0,blocksize)
!implicit none
!integer, intent(in) :: length
!integer, intent(in) :: nblocks
!integer, intent(in) :: iblock
!integer, intent(out) :: idx0
!integer, intent(out) :: blocksize
!
!integer n1,n2
!logical lerr
!
!lerr=.false.
!if (length.lt.1) lerr=.true.
!if (nblocks.lt.1) lerr=.true.
!if (iblock.lt.1.or.iblock.gt.nblocks) lerr=.true.
!if (lerr) then
!  write(*,*)
!  write(*,'("Error(idxofs) : wrong input arguments")')
!  write(*,'("  length : ",I8)')length
!  write(*,'("  nblocks : ",I8)')nblocks
!  write(*,'("  iblock : ",I8)')iblock
!  write(*,*)
!  call pstop
!endif
!n1=length/nblocks
!n2=mod(length,nblocks)
!if (n1.eq.0) then
!  if (iblock.le.n2) then
!    idx0=iblock-1
!    blocksize=1
!  else
!    idx0=-1
!    blocksize=0
!  endif
!else
!  if (iblock.le.n2) then
!    idx0=(n1+1)*(iblock-1)
!    blocksize=n1+1
!  else
!    idx0=(n1+1)*n2+n1*(iblock-n2-1)
!    blocksize=n1
!  endif
!endif
!return
!end subroutine
!
!subroutine idxloc(length,nblocks,idxg,iblock,idxl)
!implicit none
!integer, intent(in) :: length
!integer, intent(in) :: nblocks
!integer, intent(in) :: idxg
!integer, intent(out) :: iblock
!integer, intent(out) :: idxl
!
!integer n1,n2
!logical lerr
!
!lerr=.false.
!if (length.lt.1) lerr=.true.
!if (nblocks.lt.1) lerr=.true.
!if (idxg.lt.1.or.idxg.gt.length) lerr=.true.
!if (lerr) then
!  write(*,*)
!  write(*,'("Error(idxloc) : wrong input arguments")')
!  write(*,'("  length : ",I8)')length
!  write(*,'("  nblocks : ",I8)')nblocks
!  write(*,'("  idxg : ",I8)')idxg
!  write(*,*)
!  call pstop
!endif
!n1=length/nblocks
!n2=mod(length,nblocks)
!if (idxg.le.(n1+1)*n2) then
!  iblock=(idxg-1)/(n1+1)+1
!  idxl=mod(idxg-1,n1+1)+1
!else
!  iblock=n2+(idxg-(n1+1)*n2-1)/n1+1
!  idxl=mod(idxg-(n1+1)*n2-1,n1)+1
!endif
!return
!end subroutine
!
!subroutine idxglob(length,nblocks,iblock,idxl,idxg)
!implicit none
!integer, intent(in) :: length
!integer, intent(in) :: nblocks
!integer, intent(in) :: iblock
!integer, intent(in) :: idxl
!integer, intent(out) :: idxg
!
!integer n1,n2
!logical lerr
!
!lerr=.false.
!if (length.lt.1) lerr=.true.
!if (nblocks.lt.1) lerr=.true.
!if (iblock.lt.1.or.iblock.gt.nblocks) lerr=.true.
!if (idxl.lt.1) lerr=.true.
!if (lerr) then
!  write(*,*)
!  write(*,'("Error(idxglob) : wrong input arguments")')
!  write(*,'("  length : ",I8)')length
!  write(*,'("  nblocks : ",I8)')nblocks
!  write(*,'("  iblock : ",I8)')iblock
!  write(*,'("  idxl : ",I8)')idxl
!  write(*,*)
!  call pstop
!endif
!n1=length/nblocks
!n2=mod(length,nblocks)
!lerr=.false.
!if (iblock.le.n2.and.(idxl.gt.(n1+1))) lerr=.true.
!if (iblock.gt.n2.and.(idxl.gt.n1)) lerr=.true.
!if (lerr) then
!  write(*,*)
!  write(*,'("Error(idxloc) : local index out of boundary")')
!  write(*,'("  iblock : ",I8)')iblock
!  write(*,'("  idxl : ",I8)')idxl
!  write(*,*)
!  call pstop
!endif  
!if (iblock.le.n2) then
!  idxg=(iblock-1)*(n1+1)+idxl
!else
!  idxg=n2*(n1+1)+(iblock-n2-1)*n1+idxl
!endif
!return
!end subroutine
!
!  
!  
!
!
end module