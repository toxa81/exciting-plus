subroutine init_gntuju(lmaxexp)
use modmain
implicit none
integer, intent(in) :: lmaxexp
call getmaxgnt(lmaxexp,ngntujumax)
write(*,*)'x=',mpi_grid_x,'ngntujumax=',ngntujumax
call mpi_grid_barrier
if (allocated(ngntuju)) deallocate(ngntuju)
allocate(ngntuju(natmcls,ngvecme))
ngntuju=0
if (allocated(igntuju)) deallocate(igntuju)
allocate(igntuju(4,ngntujumax,natmcls,ngvecme))
igntuju=0
if (allocated(gntuju)) deallocate(gntuju)
allocate(gntuju(ngntujumax,natmcls,ngvecme))
gntuju=zzero
call gengntuju(lmaxexp)
return
end