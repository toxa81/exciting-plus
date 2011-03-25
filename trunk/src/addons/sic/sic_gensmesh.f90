subroutine sic_gensmesh
use modmain
use mod_sic
implicit none
integer itp,lm
real(8) a,tp(2)
real(8), allocatable :: glw(:),glx(:)
integer nt,np,it,ip
!
nt=lmaxwan+1
np=2*lmaxwan+1
s_ntp=nt*np
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
! spherical grid with Gauss quadrature
allocate(glw(nt),glx(nt))
call gaulegf(-1.d0,1.d0,glx,glw,nt)
itp=0
do it=1,nt
  do ip=1,np
    itp=itp+1
    s_tp(1,itp)=acos(glx(it))
    s_tp(2,itp)=twopi*(ip-1)/np
    s_tpw(itp)=twopi*glw(it)/np
    s_x(:,itp)=(/sin(s_tp(1,itp))*cos(s_tp(2,itp)),&
                 sin(s_tp(1,itp))*sin(s_tp(2,itp)),&
                 glx(it)/)
  enddo
enddo
deallocate(glx,glw)
! generate spherical harmonics
do itp=1,s_ntp
  call genrlm(lmaxwan,s_tp(1,itp),s_rlmf(1,itp))
  call genylm(lmaxwan,s_tp(1,itp),s_ylmf(1,itp))
  do lm=1,lmmaxwan
    s_rlmb(itp,lm)=s_rlmf(lm,itp)*s_tpw(itp)
    s_ylmb(itp,lm)=dconjg(s_ylmf(lm,itp))*s_tpw(itp) 
  enddo
enddo
! generate Lebedev mesh for muffin-tins
! TODO: change to Gauss
if (allocated(mt_spx)) deallocate(mt_spx)
allocate(mt_spx(3,mt_ntp))
if (allocated(mt_tpw)) deallocate(mt_tpw)
allocate(mt_tpw(mt_ntp))
if (allocated(mt_ylmf)) deallocate(mt_ylmf)
allocate(mt_ylmf(lmmaxvr,mt_ntp))
call leblaik(mt_ntp,mt_spx,mt_tpw)
do itp=1,mt_ntp
  mt_tpw(itp)=mt_tpw(itp)*fourpi
  call sphcrd(mt_spx(:,itp),a,tp)
  call genylm(lmaxvr,tp,mt_ylmf(1,itp)) 
enddo
!call nfsft_init(lmaxwan,s_ntp,stp)
return
end

