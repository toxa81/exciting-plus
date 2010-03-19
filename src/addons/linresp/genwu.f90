subroutine genwu(iwloc,vcgq,chi0,megqwan1,vscr)
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
!allocate(mtrx1(ngvecme,ngvecme))

!! rpa kernel
!mtrx1=dcmplx(0.d0,0.d0)
!do i=1,ngvecme
!  mtrx1(i,i)=vcgq(i)**2
!enddo
!! compute matrix epsilon=1-chi0*v
!epsilon=dcmplx(0.d0,0.d0)
!do i=1,ngvecme
!  epsilon(i,i)=dcmplx(1.d0,0.d0)
!enddo
!call zgemm('N','N',ngvecme,ngvecme,ngvecme,dcmplx(-1.d0,0.d0),chi0,ngvecme,mtrx1,ngvecme,&
!  dcmplx(1.d0,0.d0),epsilon,ngvecme)
!! invert epsilon matrix
!call invzge(epsilon,ngvecme)
!! compute chi=epsilon^-1 * chi0
!call zgemm('N','N',ngvecme,ngvecme,ngvecme,dcmplx(1.d0,0.d0),epsilon,ngvecme,chi0,ngvecme,&
!  dcmplx(0.d0,0.d0),mtrx1,ngvecme)

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

!! q-vector in lattice coordinates
!do i=1,3
!  vq0l(i)=1.d0*(ivq0m(i))/ngridk(i)
!enddo
!! find G-vector which brings q0 to first BZ
!vgq0l(:)=floor(vq0l(:))
!! reduce q0 vector to first BZ
!vq0rl(:)=vq0l(:)-vgq0l(:)
!! get Cartesian coordinates of q-vector and reduced q-vector
!call r3mv(bvec,vq0l,vq0c)
!call r3mv(bvec,vq0rl,vq0rc)

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
  expiqt=exp(zi*dot_product(vq0c,vtc))
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
! compute bare u
!ubare=zzero
if (mpi_grid_x(dim_k).eq.0.and.iwloc.eq.1) then
  do n1=1,nwann
    do n2=1,nwann
      do ig1=1,ngvecme
        ubarewan(n1,n2)=ubarewan(n1,n2)+&
          dconjg(megqwan(idxmegqwan(n1,n1,0,0,0),ig1))*(vcgq(ig1)**2)*&
          megqwan(idxmegqwan(n2,n2,0,0,0),ig1)
      enddo
    enddo
  enddo
endif
deallocate(epsilon,zm1)
return 
end