subroutine gensmesh_uniform
use modmain
use mod_sic
implicit none
integer itp,lm

s_ntp=lmmaxwan+500
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
call sphcover(s_ntp,s_tp)
do itp=1,s_ntp
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

