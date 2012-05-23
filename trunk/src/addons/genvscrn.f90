subroutine genvscrn(iq,chi0,krnl,vscrn,epsilon)
use modmain
use mod_addons_q
use mod_expigqr
use mod_linresp
implicit none
! arguments
integer, intent(in) :: iq
complex(8), intent(in) :: chi0(ngq(iq),ngq(iq))
complex(8), intent(in) :: krnl(ngq(iq),ngq(iq))
complex(8), intent(out) :: vscrn(ngq(iq),ngq(iq))
complex(8), intent(out) :: epsilon(ngq(iq),ngq(iq))
! local variables
!real(8), allocatable :: vcgq(:)
integer ig,ig1,ig2

call papi_timer_start(pt_vscrn)

! compute screened Coulomb potential using "symmetrized" dielectric function
!allocate(vcgq(ngq(iq)))
!do ig=1,ngq(iq)
!  vcgq(ig)=sqrt(vhgq(ig,iq))
!enddo
!do ig1=1,ngq(iq)
!  do ig2=1,ngq(iq)
!    epsilon(ig1,ig2)=-vcgq(ig1)*chi0(ig1,ig2)*vcgq(ig2)
!  enddo
!  epsilon(ig1,ig1)=zone+epsilon(ig1,ig1)
!enddo
!call invzge(epsilon,ngq(iq))
!do ig1=1,ngq(iq)
!  do ig2=1,ngq(iq)
!    vscrn(ig1,ig2)=vcgq(ig1)*epsilon(ig1,ig2)*vcgq(ig2)
!  enddo
!enddo
epsilon=zzero
do ig=1,ngq(iq)
  epsilon(ig,ig)=zone
enddo
call zgemm('N','N',ngq(iq),ngq(iq),ngq(iq),-zone,chi0,ngq(iq),&
  &krnl,ngq(iq),zone,epsilon,ngq(iq))
call invzge(epsilon,ngq(iq))
call zgemm('N','N',ngq(iq),ngq(iq),ngq(iq),zone,krnl,ngq(iq),epsilon,ngq(iq),&
  &zzero,vscrn,ngq(iq))
if (vq_gamma(iq)) then
  do ig=1,ngq(iq)
    if (igqig(ig,iq).eq.1) then
      vscrn(ig,ig)=epsilon(ig,ig)*fourpi*q0wt*nkptnr*omega/(twopi**3)
    endif
  enddo
  vscrn(:,:)=vscrn(:,:)/dble(nvq0)
endif
call papi_timer_stop(pt_vscrn)

return
end
