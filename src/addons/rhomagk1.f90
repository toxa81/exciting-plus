subroutine rhomagk1(ikloc,evecfv,evecsv)
use modmain
implicit none
! arguments
integer, intent(in) :: ikloc
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
! local variables
integer ist
integer is,ias,ik,ispn
integer ivg3(3),ifg3
integer ig1,ig2,ic,io1,io2,l1,l2,lm1,lm2,lm3
real(8) wo
complex(8) zt2(nspinor)
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)

complex(8), allocatable :: wfsvmt(:,:,:,:,:)
complex(8), allocatable :: wfsvit(:,:,:)
real(8), allocatable :: rhomtk(:,:,:,:)
complex(8), allocatable :: rhoitk(:,:)

ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)

allocate(wfsvmt(lmmaxvr,nufrmax,natmtot,nspinor,nstsv))
allocate(wfsvit(ngkmax,nspinor,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
! find the matching coefficients
do ispn=1,nspnfv
  call match(ngk(ispn,ik),gkc(:,ispn,ikloc),tpgkc(:,:,ispn,ikloc), &
   sfacgk(:,:,ispn,ikloc),apwalm(:,:,:,:,ispn))
end do
! generate wave functions in muffin-tins
call genwfsvmt(lmaxvr,lmmaxvr,ngk(1,ik),evecfv,evecsv,apwalm,wfsvmt)
! generate wave functions in interstitial
call genwfsvit(ngk(1,ik),evecfv,evecsv,wfsvit)
allocate(rhomtk(nrmtmax,lmmaxvr,natmtot,2))
rhomtk=0.d0
do ias=1,natmtot
  ic=ias2ic(ias)
  is=ias2is(ias)
  do lm3=1,lmmaxvr
    do l1=0,lmaxvr
      do io1=1,nufr(l1,is)
        do l2=0,lmaxvr
          do io2=1,nufr(l2,is)
            zt2=zzero
            do ist=1,nstsv
              wo=wkpt(ik)*occsv(ist,ik)
              if (abs(wo).gt.epsocc) then
                do lm1=l1**2+1,(l1+1)**2
                  do lm2=l2**2+1,(l2+1)**2
                    do ispn=1,nspinor
                      zt2(ispn)=zt2(ispn)+dconjg(wfsvmt(lm1,io1,ias,ispn,ist))*&
                        wfsvmt(lm2,io2,ias,ispn,ist)*sv_gntyry(lm3,lm2,lm1)*wo
                    enddo
                  enddo
                enddo
              endif
            enddo !ist
            if (spinpol) then
              rhomtk(:,lm3,ias,1)=rhomtk(:,lm3,ias,1)+dreal(zt2(1))*ufr(:,l1,io1,ic)*&
                ufr(:,l2,io2,ic)
              rhomtk(:,lm3,ias,2)=rhomtk(:,lm3,ias,2)+dreal(zt2(2))*ufr(:,l1,io1,ic)*&
                ufr(:,l2,io2,ic)
            endif
          enddo
        enddo
      enddo
    enddo
  enddo !lm3
enddo !ias
do ias=1,natmtot
  do lm3=1,lmmaxvr
    rhomt(lm3,:,ias)=rhomt(lm3,:,ias)+rhomtk(:,lm3,ias,1)+rhomtk(:,lm3,ias,2)
    magmt(lm3,:,ias,1)=magmt(lm3,:,ias,1)+rhomtk(:,lm3,ias,1)-rhomtk(:,lm3,ias,2)
  enddo
enddo
deallocate(rhomtk,apwalm,wfsvmt)

allocate(rhoitk(ngrtot,2))
rhoitk=zzero
do ig1=1,ngk(1,ik)
  do ig2=1,ngk(1,ik)
    ivg3(:)=ivg(:,igkig(ig2,1,ikloc))-ivg(:,igkig(ig1,1,ikloc))
    ifg3=igfft(ivgig(ivg3(1),ivg3(2),ivg3(3)))
    do ist=1,nstsv
      wo=wkpt(ik)*occsv(ist,ik)
      if (abs(wo).gt.epsocc) then
        do ispn=1,nspinor
          rhoitk(ifg3,ispn)=rhoitk(ifg3,ispn)+&
            dconjg(wfsvit(ig1,ispn,ist))*wfsvit(ig2,ispn,ist)*wo/omega
        enddo
      endif
    enddo
  enddo
enddo
do ispn=1,nspinor
  call zfftifc(3,ngrid,1,rhoitk(1,ispn))
enddo
rhoir(:)=rhoir(:)+dreal(rhoitk(:,1))+dreal(rhoitk(:,2))
magir(:,1)=magir(:,1)+dreal(rhoitk(:,1))-dreal(rhoitk(:,2))
deallocate(rhoitk,wfsvit)
return
end subroutine

!rho(r) = wf^{*}(r)*wf(r) = sum_{L1,L2,io1,io2} 
!  A_{L1}^{*,io1}u_{l1}^{io1}(r)*Y_{L1}^{*} *
!  A_{L2}^{io2}u_{l2}^{io2}(r)*Y_{L2} = 
!  sum_{L3} rho_{L3}(r)R_{L3}
!  
!rho_{L3}(r) = sum_{L1,L2,io1,io2} 
!  A_{L1}^{*,io1}u_{l1}^{io1}(r)*
!  A_{L2}^{io2}u_{l2}^{io2}(r) < Y_{L1} | R_{L3} | Y_{L2} > 
!  
!  do lm3=1,lmmaxvr
!    do l1=0,lmaxvr
!      do io1=1,nufr(l1,is)
!        do l2=1,lmaxvr
!          do io2=1,nufr(l2,is)
!            zt2=zzero
!            do lm1=l1**2+1,(l1+1)**2
!              do lm2=l2**2+1,(l2+1)**2
!                zt2=zt2+dconjg(A(lm1,io1))*A(lm2,io2)*gnt(lm2,lm1,lm3)
!              enddo
!            enddo
!            rho(:,lm3)=rho(:,lm3)+zt2*u(:,l1,io1)*u(:,l2,io2)
!          enddo
!        enddo
!      enddo
!    enddo
!  enddo
!            
!            
!            
!rgo(r) = wf^{*}(r)*wf(r) = 
!  sum_{G1,G2} exp^{-iG1r}C_{G1}^{*} exp^{iG2r}C_{G2} = 
!  sum_{G2} exp^{iG3r}rho_{G3}
!
!rho_{G3} = sum_{G1,G2}C_{G1}^{*}C_{G2} <G1|-G3|G2>
!  -G3+G2-G1=0 
!  G3=G2-G1   


