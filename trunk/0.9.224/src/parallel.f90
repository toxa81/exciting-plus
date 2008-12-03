subroutine barrier
#ifdef _MPI_
use mpi
integer ierr
call mpi_barrier(MPI_COMM_WORLD,ierr)
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


