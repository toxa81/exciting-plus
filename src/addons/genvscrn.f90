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
complex(8), allocatable :: krnl(:,:)
complex(8), allocatable :: chi(:,:)
real(8), allocatable :: vcgq(:)
integer ig,ig1,ig2
allocate(epsilon(ngvecme,ngvecme))
allocate(krnl(ngvecme,ngvecme))
allocate(chi(ngvecme,ngvecme))
allocate(vcgq(ngvecme))
!do ig=1,ngvecme
!  vcgq(ig)=sqrt(vhgq(ig,iq))
!enddo
! compute screened Coulomb potential using "symmetrized" dielectric function
!do ig1=1,ngvecme
!  do ig2=1,ngvecme
!    epsilon(ig1,ig2)=-vcgq(ig1)*chi0(ig1,ig2)*vcgq(ig2)
!  enddo
!  epsilon(ig1,ig1)=zone+epsilon(ig1,ig1)
!enddo
!call invzge(epsilon,ngvecme)
!do ig1=1,ngvecme
!  do ig2=1,ngvecme
!    vscrn(ig1,ig2)=vcgq(ig1)*epsilon(ig1,ig2)*vcgq(ig2)
!  enddo
!enddo
krnl=zzero
do ig=1,ngvecme
  krnl(ig,ig)=vhgq(ig,iq)
enddo
epsilon=zzero
do ig=1,ngvecme
  epsilon(ig,ig)=zone
enddo
call zgemm('N','N',ngvecme,ngvecme,ngvecme,-zone,chi0,ngvecme,&
  krnl,ngvecme,zone,epsilon,ngvecme)
call invzge(epsilon,ngvecme)
call zgemm('N','N',ngvecme,ngvecme,ngvecme,zone,epsilon,ngvecme,chi0,&
  ngvecme,zzero,chi,ngvecme)
! compute screened Coulomb potential: vscr=vbare+vbare*chi*vbare
vscrn=krnl
! compute zm2=chi*v
call zgemm('N','N',ngvecme,ngvecme,ngvecme,zone,chi,ngvecme,krnl,ngvecme,&
  zzero,epsilon,ngvecme)
! compute krnl_scr=v*zm2
call zgemm('N','N',ngvecme,ngvecme,ngvecme,zone,krnl,ngvecme,epsilon,ngvecme,&
  zone,vscrn,ngvecme)

if (all(vqm(:,iq).eq.0)) then
  do ig1=1,ngvecme
    do ig2=1,ngvecme
      if (igqig(ig1,iq).eq.1.and.ig1.eq.ig2) then
        vscrn(ig1,ig1)=vscrn(ig1,ig1)*aq0(iq)
      else
        vscrn(ig1,ig2)=vscrn(ig1,ig2)*0.125d0
      endif
    enddo !ig2      
  enddo !ig1
endif
deallocate(epsilon,vcgq)
deallocate(krnl,chi)
return
end
