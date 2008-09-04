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
   