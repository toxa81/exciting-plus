program main
use modp
implicit none
integer comm

complex(8) dat(4)
integer idx0,size,x,offs,loc,glob
integer ik

call mpi_world_initialize
call mpi_cart_initialize((/2,4/))

loc=1
glob=cart_map(17,1,loc=loc)
write(*,*)'cart_x=',cart_x,'glob=',glob

call mpi_cart_finalize
call mpi_world_finalize


return
end
