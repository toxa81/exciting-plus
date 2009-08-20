program main
use modp
implicit none
integer comm

complex(8) dat(4)


call mpi_world_initialize
call mpi_cart_initialize((/4,2/))

dat=0.d0
if (cart_root((/1/))) dat=cart_x(2)+1
call cart_bcast(dat,4,dims=(/1/))
write(*,*)'x=',cart_x,'dat=',dat

call pstop(-100)
call mpi_cart_finalize
call mpi_world_finalize


return
end
