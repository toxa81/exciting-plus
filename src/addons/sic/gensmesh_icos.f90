subroutine gensmesh_icos
use modmain
use mod_sic
implicit none
real R(0:19,3,3), v(0:11,3),vector(3)
integer resolution,n,itp
real(8) a

resolution=sic_smesh_n
n=2*resolution*(resolution-1)
n=20*n+12

if (mpi_grid_root()) then
  write(*,'("[gensmesh_icos] number of points: ",I4)')n
endif

call compute_matrices(R)
call compute_corners(v)

s_ntp=n
if (allocated(s_tp)) deallocate(s_tp)
allocate(s_tp(2,s_ntp))
if (allocated(s_x)) deallocate(s_x)
allocate(s_x(3,s_ntp))
if (allocated(s_tpw)) deallocate(s_tpw)
allocate(s_tpw(s_ntp))
do itp=1,s_ntp
  n=itp-1
  call pixel2vector(n,resolution,R,v,vector)
  s_x(:,itp)=vector(:)
  call sphcrd(s_x(1,itp),a,s_tp(1,itp))
  s_tpw(itp)=fourpi/s_ntp
enddo
return
end subroutine



