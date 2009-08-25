program main
use modp
implicit none
integer comm

complex(8) dat(4)
integer idx0,size,x,offs,loc,glob
integer ik

call mpi_world_initialize
call mpi_cart_initialize((/2,4/))

call cart_set_map_size((/17,2/))
x=4
size=cart_map(2,x=x,offs=offs)
ik=14
loc=cart_map(1,glob=ik)
write(*,*)'x=',cart_x,'offs=',offs,'size=',size,'loc=',loc,'x=',x

call mpi_cart_finalize
call mpi_world_finalize


return
end
