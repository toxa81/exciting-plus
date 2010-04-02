subroutine genuscrn(iwloc,vcgq,chi0,megqwan1,vscr)
use modmain
implicit none
integer, intent(in) :: iwloc
real(8), intent(in) :: vcgq(ngvecme)
complex(8), intent(in) :: chi0(ngvecme,ngvecme)
complex(8), intent(in) :: megqwan1(ngvecme,nwann)
complex(8), intent(out) :: vscr(ngvecme,ngvecme)

complex(8), allocatable :: epsilon(:,:)
complex(8), allocatable :: zm1(:,:)
complex(8) expiqt
real(8) vtc(3)
integer ig1,ig2,n1,n2
integer ntloc,itloc,it

allocate(epsilon(ngvecme,ngvecme))

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
    vscr(ig1,ig2)=vcgq(ig1)*epsilon(ig1,ig2)*vcgq(ig2)
  enddo
enddo
! uscrn(n,n',T,w)=<W_n,W_n|\chi(w)|W_{n',T},W_{n',T}> = 
!  \sum_{q,G,G'} megqwan^{*}_{n,G} \chi_{G,G'}(w) megqwan_{n',G'} exp^{iqT}
!  M1_{n,G'}(w)=\sum_{G} megqwan^{*}_{n,G} \chi_{G,G'}(w)
!  uscrn(n,n',T,w)=\sum_{G'} M1_{n,G'}(w) megqwan_{n',G'} exp^{iqT}
allocate(zm1(nwann,ngvecme))
ntloc=mpi_grid_map(ntr_uscrn,dim_b)
do itloc=1,ntloc
  it=mpi_grid_map(ntr_uscrn,dim_b,loc=itloc)
  vtc(:)=vtl_uscrn(1,it)*avec(:,1)+vtl_uscrn(2,it)*avec(:,2)+&
         vtl_uscrn(3,it)*avec(:,3)
  expiqt=exp(-zi*dot_product(vq0c,vtc))
  call zgemm('C','N',nwann,ngvecme,ngvecme,zone,megqwan1,ngvecme,&
    vscr,ngvecme,zzero,zm1,nwann)
  call zgemm('N','N',nwann,nwann,ngvecme,expiqt,zm1,nwann,megqwan1,ngvecme,&
    zone,uscrnwan(1,1,it,iwloc),nwann)
enddo      
! compute screened u
!do ig1=1,ngvecme
!  do ig2=1,ngvecme
!    do n1=1,nwann
!      do n2=1,nwann
!!        uscrnwan(n1,n2,ivtit_uscrn(0,0,0),iwloc)=uscrnwan(n1,n2,ivtit_uscrn(0,0,0),iwloc)+&
!!          dconjg(megqwan(idxmegqwan(n1,n1,0,0,0),ig1))*vscr(ig1,ig2)*&
!!          megqwan(idxmegqwan(n2,n2,0,0,0),ig2)
!        uscrnwan(n1,n2,ivtit_uscrn(0,0,0),iwloc)=uscrnwan(n1,n2,ivtit_uscrn(0,0,0),iwloc)+&
!          dconjg(megqwan1(ig1,n1))*vscr(ig1,ig2)*megqwan1(ig2,n2)
!      enddo
!    enddo
!  enddo
!enddo
deallocate(epsilon,zm1)
return 
end