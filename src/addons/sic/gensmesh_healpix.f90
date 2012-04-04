subroutine gensmesh_healpix
use modmain
use mod_sic
implicit none
integer nside,itp,lm
nside=12

s_ntp=12*nside**2
if (allocated(s_tp)) deallocate(s_tp)
allocate(s_tp(2,s_ntp))
if (allocated(s_x)) deallocate(s_x)
allocate(s_x(3,s_ntp))
if (allocated(s_tpw)) deallocate(s_tpw)
allocate(s_tpw(s_ntp))
if (allocated(s_rlmf)) deallocate(s_rlmf)
allocate(s_rlmf(lmmaxwan,s_ntp))
if (allocated(s_ylmf)) deallocate(s_ylmf)
allocate(s_ylmf(lmmaxwan,s_ntp))
if (allocated(s_rlmb)) deallocate(s_rlmb)
allocate(s_rlmb(s_ntp,lmmaxwan))
if (allocated(s_ylmb)) deallocate(s_ylmb)
allocate(s_ylmb(s_ntp,lmmaxwan))
do itp=1,s_ntp
  call get_healpix_tp(nside,itp-1,s_tp(1,itp),s_tp(2,itp))
  s_x(:,itp)=(/sin(s_tp(1,itp))*cos(s_tp(2,itp)),&
               sin(s_tp(1,itp))*sin(s_tp(2,itp)),&
               cos(s_tp(1,itp))/)
  s_tpw(itp)=fourpi/s_ntp
enddo

do itp=1,s_ntp
  call genrlm(lmaxwan,s_tp(1,itp),s_rlmf(1,itp))
  call genylm(lmaxwan,s_tp(1,itp),s_ylmf(1,itp))
  do lm=1,lmmaxwan
    s_rlmb(itp,lm)=s_rlmf(lm,itp)*s_tpw(itp)
    s_ylmb(itp,lm)=dconjg(s_ylmf(lm,itp))*s_tpw(itp) 
  enddo
enddo

return
end subroutine

subroutine get_healpix_tp(nside,ipix,theta,phi)
use modmain
implicit none
integer, intent(in) :: nside
integer, intent(in) :: ipix  ! in [0,npix) interval
real(8), intent(out) :: theta
real(8), intent(out) :: phi
integer npix,nl2,ncap,iring,iphi,ip,nl4
real(8) fodd

npix=12*nside**2
nl2=2*nside
ncap=nl2*(nside-1)
! north polar cap
if (ipix.lt.ncap) then
  iring=nint(sqrt((ipix+1)*0.5d0))
  iphi=ipix-2*iring*(iring-1)
  theta=acos(1.d0-(iring/dble(nside))**2/3.d0)
  phi=(dble(iphi)+0.5d0)*0.5d0*pi/iring
! equatorial region
elseif (ipix.lt.npix-ncap) then
  ip=ipix-ncap
  nl4=4*nside
  iring=int(dble(ip)/nl4)+nside
  iphi=iand(ip,nl4-1)
  fodd=0.5d0*(iand(iring+nside+1,1))
  theta=acos((nl2-iring)/(1.5d0*nside))
  phi=(dble(iphi)+fodd)*0.5d0*pi/nside
! south polar cap
else 
  ip=npix-ipix
  iring=nint(sqrt(ip*0.5d0))
  iphi=2*iring*(iring+1)-ip
  theta=acos((iring/dble(nside))**2/3.d0-1.d0)
  phi=(dble(iphi)+0.5d0)*0.5d0*pi/iring
endif
return
end subroutine


