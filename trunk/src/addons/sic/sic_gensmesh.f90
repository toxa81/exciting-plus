subroutine sic_gensmesh
use modmain
use mod_sic
implicit none
integer itp,lm
real(8) a,tp(2)
real(8), allocatable :: tp1(:,:)
!
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
if (allocated(s_spx)) deallocate(s_spx)
allocate(s_spx(3,s_ntp))
! Lebedev-Laikov mesh
call leblaik(s_ntp,s_spx,s_tpw)
! get (theta,phi) of each spx vector and generate spherical harmonics
do itp=1,s_ntp
  s_tpw(itp)=s_tpw(itp)*fourpi
  call sphcrd(s_spx(:,itp),a,tp)
  call genrlm(lmaxwan,tp,s_rlmf(1,itp))
  call genylm(lmaxwan,tp,s_ylmf(1,itp))
  do lm=1,lmmaxwan
    s_rlmb(itp,lm)=s_rlmf(lm,itp)*s_tpw(itp)
    s_ylmb(itp,lm)=dconjg(s_ylmf(lm,itp))*s_tpw(itp) 
  enddo
enddo

! uniform cover
!allocate(tp1(2,s_ntp))
!call sphcover(s_ntp,tp1)
!do itp=1,s_ntp
!  s_tpw(itp)=fourpi/s_ntp
!  s_spx(:,itp)=(/sin(tp1(1,itp))*cos(tp1(2,itp)), &
!                 sin(tp1(1,itp))*sin(tp1(2,itp)), &
!                 cos(tp1(1,itp))/)
!  call genrlm(lmaxwan,tp1(1,itp),s_rlmf(1,itp))
!  call genylm(lmaxwan,tp1(1,itp),s_ylmf(1,itp))
!  do lm=1,lmmaxwan
!    s_rlmb(itp,lm)=s_rlmf(lm,itp)*s_tpw(itp)
!    s_ylmb(itp,lm)=dconjg(s_ylmf(lm,itp))*s_tpw(itp) 
!  enddo
!enddo
!deallocate(tp1)

return
end

!
!subroutine test1(ntp,lmax)
!use modmain
!implicit none
!integer, intent(in) :: ntp
!integer, intent(in) :: lmax
!integer itp,lm,lmmax,itp1
!real(8) a,tp(2)
!complex(8) z1
!real(8), allocatable :: spx(:,:)
!real(8), allocatable :: tpw(:)
!complex(8), allocatable :: ylmf(:,:)
!complex(8), allocatable :: ylmb(:,:)
!
!lmmax=(lmax+1)**2
!
!allocate(spx(3,ntp))
!allocate(tpw(ntp))
!allocate(ylmf(lmmax,ntp))
!allocate(ylmb(ntp,lmmax))
!! Lebedev-Laikov mesh
!call leblaik(ntp,spx,tpw)
!! get (theta,phi) of each spx vector and generate spherical harmonics
!do itp=1,ntp                   
!  tpw(itp)=tpw(itp)*fourpi
!  call sphcrd(spx(1,itp),a,tp)
!  call genylm(lmax,tp,ylmf(1,itp))
!  do lm=1,lmmax
!    ylmb(itp,lm)=dconjg(ylmf(lm,itp))*tpw(itp) 
!  enddo  
!enddo 
!a=0.d0
!do itp=1,ntp
!  do itp1=1,ntp
!    z1=zzero
!    do lm=1,lmmax
!      z1=z1+ylmb(itp,lm)*ylmf(lm,itp1)
!    enddo
!    !if (itp.eq.itp1) z1=z1-zone
!    !a=max(a,abs(z1))
!    if (itp.ne.itp1) a=max(a,abs(z1))
!  enddo
!enddo
!if (mpi_grid_root()) then
!  write(*,'("[test1] l = ",I4," ntp = ",I4)')lmax,ntp
!  write(*,'("[test1] Lebedev quadrature completeness error : ",G18.10)')a
!  write(*,*)
!endif
!deallocate(spx,tpw,ylmf,ylmb)
!return
!end

