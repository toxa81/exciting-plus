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
