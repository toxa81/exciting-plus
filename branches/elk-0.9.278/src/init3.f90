subroutine init3
use modmain
implicit none
integer, allocatable :: dims(:)

if (task.eq.0.or.task.eq.1) then
  allocate(dims(2))
  if (nproc.le.nkpt) then
    dims(dim_nkpt)=nproc
    dims(dim_xtra)=1
  else
    dims(dim_nkpt)=nkpt
    dims(dim_xtra)=nproc/nkpt
  endif
  call mpi_cart_initialize(dims)
  nkptloc=cart_map(nkpt,dim_nkpt)
  deallocate(dims)
endif
return
end