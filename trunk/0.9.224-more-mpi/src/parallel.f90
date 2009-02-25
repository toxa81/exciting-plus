subroutine pinit_cart
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
#ifdef _MPI_
integer i,ierr
logical, allocatable :: mpi_periods(:)
!if (.not.lmpi_init) then
!  call mpi_init(ierr)
!  comm_world=MPI_COMM_WORLD
!  call mpi_comm_size(comm_world,nproc,ierr)
!  call mpi_comm_rank(comm_world,iproc,ierr)
!  lmpi_init=.true.
!endif
#endif
mpi_ndims=1
if (task.eq.400.or.task.eq.401.or.task.eq.402) then
  mpi_ndims=3
endif
if (allocated(mpi_dims)) deallocate(mpi_dims)
allocate(mpi_dims(mpi_ndims))
if (allocated(mpi_x)) deallocate(mpi_x)
allocate(mpi_x(mpi_ndims))
#ifdef _MPI_
mpi_dims=nproc
if (task.eq.400.or.task.eq.401) then
  if (nproc.le.nkptnr) then
    mpi_dims=(/nproc,1,1/)
  else
    i=nproc/nkptnr
    if (i.le.nvq0) then
      mpi_dims=(/nkptnr,1,i/)
    else
      mpi_dims=(/nkptnr,nproc/(nkptnr*nvq0),nvq0/)
    endif
  endif
endif
if (task.eq.402) then
  if (nproc.le.nvq0) then
    mpi_dims=(/1,1,nproc/)
  else
    i=nproc/nvq0
    if (i.le.nfxca) then
      mpi_dims=(/1,i,nvq0/)
    else
      mpi_dims=(/nproc/(nvq0*nfxca),nfxca,nvq0/)
    endif
  endif
endif
allocate(mpi_periods(mpi_ndims))
mpi_periods=.false.
call mpi_cart_create(comm_world,mpi_ndims,mpi_dims,mpi_periods,   &
  .false.,comm_cart,ierr)
call mpi_cart_get(comm_cart,mpi_ndims,mpi_dims,mpi_periods,mpi_x, &
  ierr)
deallocate(mpi_periods)
if (mpi_ndims.eq.3) then
  call mpi_cart_sub(comm_cart,(/.true.,.false.,.false./),comm_cart_100,ierr)
  call mpi_cart_sub(comm_cart,(/.false.,.true.,.false./),comm_cart_010,ierr)
  call mpi_cart_sub(comm_cart,(/.false.,.false.,.true./),comm_cart_001,ierr)
  call mpi_cart_sub(comm_cart,(/.true.,.false.,.true./),comm_cart_101,ierr)
  call mpi_cart_sub(comm_cart,(/.false.,.true.,.true./),comm_cart_011,ierr)
  call mpi_cart_sub(comm_cart,(/.true.,.true.,.false./),comm_cart_110,ierr)
endif
#else
mpi_dims=1
mpi_x=0
#endif
return
end










subroutine barrier(comm)
#ifdef _MPI_
use mpi
#endif
implicit none
integer, intent(in) :: comm
#ifdef _MPI_
integer ierr
call mpi_barrier(comm,ierr)
#endif
return
end

subroutine pstop
#ifdef _MPI_
use mpi
integer ierr
call mpi_abort(MPI_COMM_WORLD,-1,ierr)
call mpi_finalize(ierr)
#else
stop
#endif
return
end

subroutine splitk(nkpt,nproc,nkptloc,ikptloc)
implicit none
! arguments
integer, intent(in) :: nkpt
integer, intent(in) :: nproc
integer, intent(out) :: nkptloc(0:nproc-1)
integer, intent(out) :: ikptloc(0:nproc-1,2) 
! local variables
integer i,n1,n2

! minimum number of k-points for each proc.
n1=nkpt/nproc
! remaining number of k-points which will be distributed among first n2 procs
n2=nkpt-n1*nproc
! each proc gets n1 k-points
nkptloc(:)=n1
! additionally, first n2 procs get extra point
do i=0,n2-1
  nkptloc(i)=nkptloc(i)+1
enddo 
! build index of first and last k-point for each proc
ikptloc(0,1)=1
ikptloc(0,2)=nkptloc(0)
do i=1,nproc-1
  ikptloc(i,1)=ikptloc(i-1,2)+1
  ikptloc(i,2)=ikptloc(i,1)+nkptloc(i)-1
enddo
    
return
end

integer function ikglob(ikloc)
use modmain
implicit none
integer, intent(in) :: ikloc
ikglob=ikptloc(iproc,1)+ikloc-1
return
end

integer function ikloc(ikglob)
use modmain
implicit none
integer, intent(in) :: ikglob
ikloc=ikglob-ikptloc(iproc,1)+1
return
end

integer function iknrglob(ikloc)
use modmain
implicit none
integer, intent(in) :: ikloc
iknrglob=ikptnrloc(iproc,1)+ikloc-1
return
end

integer function iknrglob2(ikloc,ipos_k)
use modmain
implicit none
integer, intent(in) :: ikloc
integer, intent(in) :: ipos_k
iknrglob2=ikptnrloc(ipos_k,1)+ikloc-1
return
end

integer function iknrloc(ikglob)
use modmain
implicit none
integer, intent(in) :: ikglob
iknrloc=ikglob-ikptnrloc(iproc,1)+1
return
end

subroutine dsync(var,n,doreduce,dobcast)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
logical, intent(in) :: doreduce
logical, intent(in) :: dobcast
integer, intent(in) :: n
real(8), intent(inout) :: var(n)

#ifdef _MPI_
integer ierr
real(8), allocatable :: tmp(:)

if (doreduce.and.iproc.eq.0) allocate(tmp(n))
if (doreduce) then
  call mpi_reduce(var,tmp,n,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
    MPI_COMM_WORLD,ierr)
  if (iproc.eq.0) var=tmp
endif
if (dobcast) call mpi_bcast(var,n,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
if (doreduce.and.iproc.eq.0) deallocate(tmp)
#endif

return
end


subroutine rsync(var,n,doreduce,dobcast)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
logical, intent(in) :: doreduce
logical, intent(in) :: dobcast
integer, intent(in) :: n
real(4), intent(inout) :: var(n)

#ifdef _MPI_
integer ierr
real(4), allocatable :: tmp(:)

if (doreduce.and.iproc.eq.0) allocate(tmp(n))
if (doreduce) then
  call mpi_reduce(var,tmp,n,MPI_REAL,MPI_SUM,0, &
    MPI_COMM_WORLD,ierr)
  if (iproc.eq.0) var=tmp
endif
if (dobcast) call mpi_bcast(var,n,MPI_REAL,0,MPI_COMM_WORLD,ierr)
if (doreduce.and.iproc.eq.0) deallocate(tmp)
#endif

return
end

subroutine zsync(var,n,doreduce,dobcast)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
logical, intent(in) :: doreduce
logical, intent(in) :: dobcast
integer, intent(in) :: n
complex(8), intent(inout) :: var(n)

#ifdef _MPI_
integer ierr
complex(8), allocatable :: tmp(:)

if (doreduce.and.iproc.eq.0) allocate(tmp(n))
if (doreduce) then
  call mpi_reduce(var,tmp,n,MPI_DOUBLE_COMPLEX,MPI_SUM,0, &
    MPI_COMM_WORLD,ierr)
  if (iproc.eq.0) var=tmp
endif
if (dobcast) call mpi_bcast(var,n,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
if (doreduce.and.iproc.eq.0) deallocate(tmp)
#endif

return
end

subroutine lsync(var,n,dobcast)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
logical, intent(in) :: dobcast
integer, intent(in) :: n
logical, intent(inout) :: var(n)

#ifdef _MPI_
integer ierr
if (dobcast) call mpi_bcast(var,n,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
#endif

return
end

subroutine isync(var,n,doreduce,dobcast)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
logical, intent(in) :: doreduce
logical, intent(in) :: dobcast
integer, intent(in) :: n
integer, intent(inout) :: var(n)

#ifdef _MPI_
integer ierr
integer, allocatable :: tmp(:)

if (doreduce.and.iproc.eq.0) allocate(tmp(n))
if (doreduce) then
  call mpi_reduce(var,tmp,n,MPI_INTEGER,MPI_SUM,0, &
    MPI_COMM_WORLD,ierr)
  if (iproc.eq.0) var=tmp
endif
if (dobcast) call mpi_bcast(var,n,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (doreduce.and.iproc.eq.0) deallocate(tmp)
#endif

return
end

logical function root_cart(idims)
use modmain
implicit none
integer, intent(in) :: idims(mpi_ndims)
logical l1
integer i
l1=.true.
do i=1,mpi_ndims
  if (idims(i).eq.1.and.mpi_x(i).ne.0) l1=.false.
enddo
root_cart=l1
return
end


subroutine i_bcast_cart(comm,val,n)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
integer, intent(in) :: comm
integer, intent(in) :: n
integer, intent(inout) :: val(n)
integer root,ierr,comm_x(mpi_ndims)
#ifdef _MPI_
comm_x=0
call mpi_cart_rank(comm,comm_x,root,ierr)
call mpi_bcast(val,n,MPI_INTEGER,root,comm,ierr)
#endif
return
end







subroutine d_reduce_cart(comm,all,val,n)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
integer, intent(in) :: comm
logical, intent(in) :: all
integer, intent(in) :: n
real(8), intent(inout) :: val(n)
real(8), allocatable :: tmp(:)
integer comm_dims(mpi_ndims),comm_x(mpi_ndims)
logical comm_periods(mpi_ndims)
integer root,rank,ierr
#ifdef _MPI_
comm_dims=-1
comm_x=-1
call mpi_cart_get(comm,mpi_ndims,comm_dims,comm_periods,comm_x,ierr)
call mpi_cart_rank(comm,comm_x,rank,ierr)
comm_x=0
call mpi_cart_rank(comm,comm_x,root,ierr)
if (all) then
  allocate(tmp(n))
  call mpi_allreduce(val,tmp,n,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  val=tmp
  deallocate(tmp)
else
  if (rank.eq.root) allocate(tmp(n))
  call mpi_reduce(val,tmp,n,MPI_DOUBLE_PRECISION,MPI_SUM,root,comm,ierr)
  if (rank.eq.root) then
    val=tmp
    deallocate(tmp)
  endif
endif
#endif
return
end

subroutine d_bcast_cart(comm,val,n)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
integer, intent(in) :: comm
integer, intent(in) :: n
real(8), intent(inout) :: val(n)
integer root,ierr,comm_x(mpi_ndims)
#ifdef _MPI_
comm_x=0
call mpi_cart_rank(comm,comm_x,root,ierr)
call mpi_bcast(val,n,MPI_DOUBLE_PRECISION,root,comm,ierr)
#endif
return
end


subroutine idxbos(length,nblocks,iblock,idx0,blocksize)
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
  write(*,'("Error(idxbos) : wrong input arguments")')
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
end

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
end

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
end



