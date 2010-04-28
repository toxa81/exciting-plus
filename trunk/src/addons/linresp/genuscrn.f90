subroutine genuscrn(iwloc,iq,ivq0m,vcgq,chi0,megqwan1,vscr)
use modmain
implicit none
integer, intent(in) :: iwloc
integer, intent(in) :: iq
integer, intent(in) :: ivq0m(3)
real(8), intent(in) :: vcgq(ngvecme)
complex(8), intent(in) :: chi0(ngvecme,ngvecme)
complex(8), intent(in) :: megqwan1(ngvecme,nwann)
complex(8), intent(out) :: vscr(ngvecme,ngvecme)

complex(8), allocatable :: epsilon(:,:)
complex(8), allocatable :: zm1(:,:)
complex(8) expiqt
real(8) vtc(3)
integer ig1,ig2
integer ntloc,itloc,it

allocate(epsilon(ngvecme,ngvecme))

if (crpa_scrn.eq.0) then
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
endif
if (crpa_scrn.eq.1.or.crpa_scrn.eq.2) then
! compute matrix 1-chi0*v
  do ig1=1,ngvecme
    do ig2=1,ngvecme
      epsilon(ig1,ig2)=-chi0(ig1,ig2)*(vcgq(ig2)**2)
    enddo
    epsilon(ig1,ig1)=epsilon(ig1,ig1)+zone
  enddo
  call invzge(epsilon,ngvecme)
  allocate(zm1(ngvecme,ngvecme))
! compute chi=epsilon^-1 * chi0
  call zgemm('N','N',ngvecme,ngvecme,ngvecme,zone,epsilon,ngvecme,chi0,&
    ngvecme,zzero,zm1,ngvecme)
! compute screened Coulomb potential: vscr=vbare+vbare*chi*vbare
  do ig1=1,ngvecme
    do ig2=1,ngvecme
      vscr(ig1,ig2)=(vcgq(ig1)**2)*zm1(ig1,ig2)*(vcgq(ig2)**2)
    enddo
  enddo
  if (crpa_scrn.eq.1) then
    do ig1=1,ngvecme
      vscr(ig1,ig1)=vscr(ig1,ig1)+(vcgq(ig1)**2)
    enddo
  endif
  deallocate(zm1)
endif
do ig1=1,ngvecme
  if (ivq0m(1).eq.0.and.ivq0m(2).eq.0.and.ivq0m(3).eq.0) then
    if (ig1.eq.1) then
      vscr(ig1,ig1)=vscr(ig1,ig1)*a0gamma(iq)
    else
      vscr(ig1,ig1)=vscr(ig1,ig1)*0.125d0
    endif      
  endif
enddo

! uscrn(n,n',T,w)=<W_n,W_n|\chi(w)|W_{n',T},W_{n',T}> = 
!  \sum_{q,G,G'} megqwan^{*}_{n,G} \chi_{G,G'}(w) megqwan_{n',G'} exp^{-iqT}
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
deallocate(epsilon,zm1)
return 
end