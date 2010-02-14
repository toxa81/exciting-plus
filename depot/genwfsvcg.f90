subroutine genwfsvcg(ngp,igpig,gp,ylmgp,sfacgp,wfsvmt,wfsvit,wfsvcg)
use modmain
implicit none
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: gp(ngkmax)
complex(8), intent(in) :: ylmgp(lmmaxvr,ngkmax)
complex(8), intent(in) :: sfacgp(ngkmax,natmtot)
complex(8), intent(in) :: wfsvmt(lmmaxvr,nrfmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfsvit(ngkmax,nspinor,nstsv)
complex(8), intent(out) :: wfsvcg(ngkmax,nspinor,nstsv)

integer ik,ig,is,ias,ir,l1,io1,lm1,j,ispn,ig1,j1
integer ivg1(3)
real(8) fr(nrmtmax),gr(nrmtmax),cf(3,nrmtmax)
real(8), allocatable :: jl(:,:)
real(8), allocatable :: uj(:,:)
complex(8) zt1

allocate(jl(nrmtmax,0:lmaxvr))
allocate(uj(0:lmaxvr,nrfmax))

wfsvcg=zzero
do ig=1,ngp
! muffin-tin part
  do ias=1,natmtot
    uj=0.d0
    is=ias2is(ias)
! generate Bessel functions j_l(|G+k|x)
    do ir=1,nrmt(is)
      call sbessel(lmaxvr,gp(ig)*spr(ir,is),jl(ir,:))
    enddo
    do l1=0,lmaxvr
      do io1=1,nrfmax
        do ir=1,nrmt(is)
          fr(ir)=urf(ir,l1,io1,ias)*jl(ir,l1)*(spr(ir,is)**2)
        enddo !ir
        call fderiv(-1,nrmt(is),spr(1,is),fr,gr,cf)
        uj(l1,io1)=gr(nrmt(is))
      enddo !io1
    enddo !l1
    do j=1,nstsv
      do ispn=1,nspinor
        zt1=zzero
        do lm1=1,lmmaxvr
          l1=lm2l(lm1)
          do io1=1,nrfmax
            zt1=zt1+dconjg(zi**l1)*ylmgp(lm1,ig)*&
              wfsvmt(lm1,io1,ias,ispn,j)*uj(l1,io1)
          enddo
        enddo
        zt1=zt1*fourpi*dconjg(sfacgp(ig,ias)) !/omega
        wfsvcg(ig,ispn,j)=wfsvcg(ig,ispn,j)+zt1
      enddo !ispn
    enddo !j
  enddo !ias
! interstitial part
  do j=1,nstsv
    do ispn=1,nspinor
      zt1=zzero
! compute \sum_{G1} u_1^{*}(G1)*\theta(G1+G)
      do ig1=1,ngp
! ivg1(:)=G1+G             
        ivg1(:)=-ivg(:,igpig(ig1))+ivg(:,igpig(ig))
!                if (ivg1(1).lt.intgv(1,1).or.ivg1(1).gt.intgv(1,2).or.&
!                    ivg1(2).lt.intgv(2,1).or.ivg1(2).gt.intgv(2,2).or.&
!                    ivg1(3).lt.intgv(3,1).or.ivg1(3).gt.intgv(3,2)) then
!                  write(*,*)
!                  write(*,'("Error(megqblhit): G-vector is outside of boundaries")')
!                  write(*,'("  G1+G : ",3I5)')ivg1
!                  write(*,'("  boundaries : ",2I5,",",2I5,",",2I5)')&
!                    intgv(1,:),intgv(2,:),intgv(3,:)
!                  write(*,*)
!                  call pstop
!                endif
                zt1=zt1+wfsvit(ig1,ispn,j)* &
                        cfunig(ivgig(ivg1(1),ivg1(2),ivg1(3)))
  
              enddo !ig1
                   !wfsvcg(ig,ispn,j)=wfsvcg(ig,ispn,j)+zt1  
             enddo
             enddo
 
enddo !ig

do j=1,nstsv
  do j1=1,nstsv
  zt1=zzero
  do ispn=1,nspinor
  do ig=1,ngp
    zt1=zt1+dconjg(wfsvcg(ig,ispn,j))*wfsvcg(ig,ispn,j1)
  enddo
  enddo
  write(*,*)'j=',j,'j1=',j1,'zt1=',zt1
enddo
enddo

!call wfprodk(ngp,igpig,wfsvmt,wfsvit,zt1)

return
end


