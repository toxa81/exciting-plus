subroutine gensmesh_lebedev
use modmain
use mod_sic
implicit none
integer itp
real(8) a

s_ntp=sic_smesh_n

if (allocated(s_tp)) deallocate(s_tp)
allocate(s_tp(2,s_ntp))
if (allocated(s_x)) deallocate(s_x)
allocate(s_x(3,s_ntp))
if (allocated(s_tpw)) deallocate(s_tpw)
allocate(s_tpw(s_ntp))
call leblaik(s_ntp,s_x,s_tpw)
do itp=1,s_ntp
  call sphcrd(s_x(1,itp),a,s_tp(1,itp))
  s_tpw(itp)=fourpi*s_tpw(itp)
enddo
return
end subroutine

