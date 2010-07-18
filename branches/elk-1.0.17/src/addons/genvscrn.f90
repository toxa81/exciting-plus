subroutine genvscrn(iq,chi0,vscrn)
use modmain
use mod_addons_q
implicit none
! arguments
integer, intent(in) :: iq
complex(8), intent(in) :: chi0(ngvecme,ngvecme)
complex(8), intent(out) :: vscrn(ngvecme,ngvecme)
! local variables
complex(8), allocatable :: epsilon(:,:)
real(8), allocatable :: vcgq(:)
integer ig,ig1,ig2
allocate(epsilon(ngvecme,ngvecme))
allocate(vcgq(ngvecme))
do ig=1,ngvecme
  vcgq(ig)=sqrt(vhgq(ig,iq))
enddo
! compute screened Coulomb potential using "symmetrized" dielectric function
do ig1=1,ngvecme
  do ig2=1,ngvecme
    epsilon(ig1,ig2)=-vcgq(ig1)*chi0(ig1,ig2)*vcgq(ig2)
  enddo
  epsilon(ig1,ig1)=zone+epsilon(ig1,ig1)
enddo
call invzge(epsilon,ngvecme)
do ig1=1,ngvecme
  do ig2=1,ngvecme
    vscrn(ig1,ig2)=vcgq(ig1)*epsilon(ig1,ig2)*vcgq(ig2)
  enddo
enddo
if (all(vqm(:,iq).eq.0)) then
  do ig=1,ngvecme
    if (igqig(ig,iq).eq.1) then
      vscrn(ig,ig)=vscrn(ig,ig)*aq0(iq)
    else
      vscrn(ig,ig)=vscrn(ig,ig)*0.125d0
    endif      
  enddo !ig
endif
deallocate(epsilon,vcgq)
return
end
