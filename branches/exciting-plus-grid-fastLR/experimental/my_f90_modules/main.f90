program main
use mod_mpi_grid
implicit none
integer comm

real(8) dat(4),tt
integer idx0,size,x,offs,loc,glob
integer ik
integer idat(4),iidat(1)

call mpi_world_initialize
call mpi_grid_initialize((/4,2/),.true.)

call pstop
!dat=0.d0
!dat(cart_x(1)+1)=1.0*cart_x(1)+1.0
!write(*,*)'cart_x=',cart_x,'dat=',dat
!call cart_barrier()

!call cart_reduce(dat(1),4,dims=(/1/),all=.true.,side=.true.)
!write(*,*)'cart_x=',cart_x,'idat=',dat
!call cart_barrier()

call mpi_grid_finalize
call mpi_world_finalize


return
end
