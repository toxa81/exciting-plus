subroutine barrier
#ifdef _MPI_
use mpi
integer ierr
call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif
return
end

subroutine barrier2(mpi_comm)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
integer, intent(in) :: mpi_comm
#ifdef _MPI_
integer ierr
call mpi_barrier(mpi_comm,ierr)
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

subroutine dsync2(idims,var,n,doreduce,dobcast)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
integer, intent(in) :: idims(mpi_ndims)
logical, intent(in) :: doreduce
logical, intent(in) :: dobcast
integer, intent(in) :: n
real(8), intent(inout) :: var(n)

#ifdef _MPI_
integer mpi_comm_tmp,i,ierr
real(8), allocatable :: tmp(:)
logical ldims(mpi_ndims)
logical ldims2(mpi_ndims)
integer rank
logical, external :: in_set

do i=1,mpi_ndims
  if (idims(i).eq.0) then
    ldims(i)=.false.
    ldims2(i)=.true.
  else 
    ldims(i)=.true.
    ldims2(i)=.false.
  endif
enddo
call mpi_cart_sub(mpi_comm_cart,ldims,mpi_comm_tmp,ierr)

if (in_set(idims)) then
  call mpi_cart_rank(mpi_comm_tmp,mpi_x,rank,ierr)
  if (doreduce.and.rank.eq.0) allocate(tmp(n))
  if (doreduce) then
    call mpi_reduce(var,tmp,n,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
      mpi_comm_tmp,ierr)
    if (rank.eq.0) var=tmp
  endif
  if (dobcast) call mpi_bcast(var,n,MPI_DOUBLE_PRECISION,0,mpi_comm_tmp,ierr)
  if (doreduce.and.rank.eq.0) deallocate(tmp)
endif

call mpi_cart_sub(mpi_comm_cart,ldims2,mpi_comm_tmp,ierr)
call mpi_bcast(var,n,MPI_DOUBLE_PRECISION,0,mpi_comm_tmp,ierr)

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

logical function in_set(idims)
use modmain
implicit none
integer, intent(in) :: idims(mpi_ndims)
logical l1
integer i
l1=.true.
do i=1,mpi_ndims
  if (idims(i).eq.0.and.mpi_x(i).ne.0) l1=.false.
enddo
in_set=l1
return
end

subroutine split_idx(l,n,j,l0,l1)
implicit none
! arguments
integer, intent(in) :: l
integer, intent(in) :: n
integer, intent(in) :: j
integer, intent(out) :: l0
integer, intent(out) :: l1 
! local variables
integer i,n1,n2
integer tmp(n),tmp2(n,2)

! minimum number of points for each segment
n1=l/n
! remaining number of points which will be distributed among first n2 segments
n2=l-n1*n
! each segment gets n1 points
tmp(:)=n1
! additionally, first n2 segments get extra point
do i=1,n2
  tmp(i)=tmp(i)+1
enddo 
! build index of first and last point for each segment
tmp2(1,1)=1
tmp2(1,2)=tmp(1)
do i=2,n
  tmp2(i,1)=tmp2(i-1,2)+1
  tmp2(i,2)=tmp2(i,1)+tmp(i)-1
enddo
l0=tmp2(j,1)
l1=tmp2(j,2)
   
return
end

subroutine d_allreduce(mpi_comm,val,n)
use mpi
implicit none
integer, intent(in) :: mpi_comm
integer, intent(in) :: n
real(8), intent(inout) :: val(n)
real(8), allocatable :: tmp(:)
integer ierr
allocate(tmp(n))
call mpi_allreduce(val,tmp,n,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm,ierr)
val=tmp
deallocate(tmp)
return
end

subroutine d_bcast(idims,val,n)
use modmain
use mpi
implicit none
integer, intent(in) :: idims(mpi_ndims)
integer, intent(in) :: n
real(8), intent(inout) :: val(n)
logical ldims(mpi_ndims)
integer i,root,mpi_comm,ierr,nsub_coord
integer, allocatable :: sub_coord(:)
nsub_coord=0
do i=1,mpi_ndims
  if (idims(i).eq.0) then
    ldims(i)=.false.
  else
    ldims(i)=.true.
    nsub_coord=nsub_coord+1
  endif
enddo
call mpi_cart_sub(mpi_comm_cart,ldims,mpi_comm,ierr)
allocate(sub_coord(nsub_coord))
sub_coord=0
call mpi_cart_rank(mpi_comm,sub_coord,root,ierr)
call mpi_bcast(val,n,MPI_DOUBLE_PRECISION,root,mpi_comm,ierr)
deallocate(sub_coord)
return
end

subroutine i_bcast(idims,val,n)
use modmain
use mpi
implicit none
integer, intent(in) :: idims(mpi_ndims)
integer, intent(in) :: n
integer, intent(inout) :: val(n)
logical ldims(mpi_ndims)
integer i,root,mpi_comm,ierr,nsub_coord
integer, allocatable :: sub_coord(:)
nsub_coord=0
do i=1,mpi_ndims
  if (idims(i).eq.0) then
    ldims(i)=.false.
  else
    ldims(i)=.true.
    nsub_coord=nsub_coord+1
  endif
enddo
call mpi_cart_sub(mpi_comm_cart,ldims,mpi_comm,ierr)
allocate(sub_coord(nsub_coord))
sub_coord=0
call mpi_cart_rank(mpi_comm,sub_coord,root,ierr)
call mpi_bcast(val,n,MPI_INTEGER,root,mpi_comm,ierr)
deallocate(sub_coord)
return
end

subroutine d_reduce(idims,all,val,n)
use modmain
use mpi
implicit none
integer, intent(in) :: idims(mpi_ndims)
logical, intent(in) :: all
integer, intent(in) :: n
real(8), intent(inout) :: val(n)
real(8), allocatable :: tmp(:)
logical ldims(mpi_ndims)
integer i,mpi_comm,ierr

do i=1,mpi_ndims
  if (idims(i).eq.0) then
    ldims(i)=.false.
  else
    ldims(i)=.true.
  endif
enddo

call mpi_cart_sub(mpi_comm_cart,ldims,mpi_comm,ierr)
if (all) then
  allocate(tmp(n))
  call mpi_allreduce(val,tmp,n,MPI_DOUBLE_PRECISION,MPI_SUM,mpi_comm,ierr)
  val=tmp
  deallocate(tmp)
else
  write(*,*)'reduce not implemented'
endif
return
end



