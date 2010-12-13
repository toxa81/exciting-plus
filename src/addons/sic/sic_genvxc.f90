subroutine sic_genvxc(vhxcmt,vhxcir,ene)
use modmain
use mod_sic
use modxcifc
implicit none
real(8), intent(out) :: vhxcmt(lmmaxvr,nmtloc,sic_orbitals%ntr,nspinor,sic_wantran%nwan)
real(8), intent(out) :: vhxcir(ngrloc,sic_orbitals%ntr,nspinor,sic_wantran%nwan)
complex(8), intent(out) :: ene(4,sic_wantran%nwan)
integer ntp,itp,lm,n,ispn,it,j
real(8), allocatable :: tp(:,:)
complex(8), allocatable :: ylm(:,:)
complex(8), allocatable :: rlmz(:,:)
real(8) rlm(lmmaxvr)
real(8), allocatable :: rlmb(:,:)
real(8), allocatable :: excwanmt(:,:,:)
real(8), allocatable :: excwanir(:,:)
complex(8), allocatable :: wfmt(:,:)
real(8), allocatable :: wfmt2(:,:,:)
real(8), allocatable :: wfir2(:,:)
real(8), allocatable :: exmt_(:,:)
real(8), allocatable :: exir_(:)
real(8), allocatable :: ecmt_(:,:)
real(8), allocatable :: ecir_(:)
real(8), allocatable :: vxmt_(:,:,:)
real(8), allocatable :: vxir_(:,:)
real(8), allocatable :: vcmt_(:,:,:)
real(8), allocatable :: vcir_(:,:)
real(8), allocatable :: wtp(:)
real(8), allocatable :: spx(:,:)
real(8) t1

! XC energy density of Wannier function
allocate(excwanmt(lmmaxvr,nmtloc,sic_orbitals%ntr))
allocate(excwanir(ngrloc,sic_orbitals%ntr))

!nt=2*lmaxvr+1
!np=2*lmaxvr+1
!ntp=nt*np
!
!allocate(w(nt))
!allocate(tp(2,ntp),wtp(ntp))
!allocate(ylm(lmmaxvr,ntp))
!allocate(rlm1(lmmaxvr,ntp))
!allocate(rlmc(ntp,lmmaxvr))
!
!itp=0
!do it=0,2*lmaxvr
!  theta=(pi/(2*lmaxvr+1))*(0.5d0+it)
!  w(it+1)=1.d0
!  do l=1,lmaxvr
!    w(it+1)=w(it+1)+2.d0*cos(2.d0*l*theta)/(1.d0-4.d0*l**2)
!  enddo
!  w(it+1)=w(it+1)*2/pi
!  do ip=0,2*lmaxvr
!    phi=(2*pi/(2*lmaxvr+1))*ip
!    itp=itp+1
!    tp(1,itp)=theta
!    tp(2,itp)=phi
!    wtp(itp)=w(it+1)*pi*pi*2/(2*lmaxvr+1)**2
!  enddo
!enddo
!do itp=1,ntp 
!  call genylm(lmaxvr,tp(1,itp),ylm(1,itp))
!  call genrlm(lmaxvr,tp(1,itp),rlm)  
!  do lm=1,lmmaxvr
!    rlmc(itp,lm)=rlm(lm)*wtp(itp)
!    rlm1(lm,itp)=zone*rlm(lm)
!  enddo
!enddo
!deallocate(w,wtp)


!nt=lmaxvr+1
!np=2*nt
!ntp=nt*np
!allocate(x(nt),w(nt))
!allocate(tp(2,ntp),wtp(ntp))
!allocate(ylm(lmmaxvr,ntp))
!call gaulegf(-1.d0,1.d0,x,w,nt)
!itp=0
!do it=1,nt
!  do ip=1,np
!    itp=itp+1
!    tp(1,itp)=acos(x(it))
!    tp(2,itp)=2*pi*(ip-1)/np
!    wtp(itp)=pi*w(it)/nt
!  enddo
!enddo
!allocate(rlmc(ntp,lmmaxvr))
!do itp=1,ntp 
!  call genylm(lmaxvr,tp(1,itp),ylm(1,itp))
!  call genrlm(lmaxvr,tp(1,itp),rlm)  
!  do lm=1,lmmaxvr
!    rlmc(itp,lm)=rlm(lm)*wtp(itp)*sin(tp(1,itp))
!  enddo
!enddo
!deallocate(x,w,wtp)

!ntp=1000
!allocate(tp(2,ntp))
!allocate(ylm(lmmaxvr,ntp))
!allocate(rlmz(lmmaxvr,ntp))
!allocate(rlmb(ntp,lmmaxvr))
!call sphcover(ntp,tp)
!do itp=1,ntp 
!  call genylm(lmaxvr,tp(1,itp),ylm(1,itp))
!  call genrlm(lmaxvr,tp(1,itp),rlm)  
!  do lm=1,lmmaxvr
!    rlmc(itp,lm)=rlm(lm)*fourpi/ntp
!    rlmz(lm,itp)=zone*rlm(lm)
!  enddo
!enddo

ntp=266
! allocate common arrays
allocate(tp(2,ntp))
allocate(ylm(lmmaxvr,ntp))
allocate(rlmz(lmmaxvr,ntp))
allocate(rlmb(ntp,lmmaxvr))
allocate(wtp(ntp))

! uniform cover
!call sphcover(ntp,tp)
!do itp=1,ntp
!  wtp(itp)=fourpi/ntp
!enddo

! Lebedev mesh
allocate(spx(3,ntp))
call leblaik(ntp,spx,wtp)                                                                   
do itp=1,ntp                   
  wtp(itp)=wtp(itp)*fourpi
  call sphcrd(spx(:,itp),t1,tp(:, itp))                                                     
enddo 
deallocate(spx)

! generate spherical harmonics
do itp=1,ntp 
  call genylm(lmaxvr,tp(1,itp),ylm(1,itp))
  call genrlm(lmaxvr,tp(1,itp),rlm)  
  do lm=1,lmmaxvr
    rlmb(itp,lm)=rlm(lm)*wtp(itp)
    rlmz(lm,itp)=zone*rlm(lm)
  enddo
enddo
deallocate(wtp,tp)

allocate(wfmt(ntp,nmtloc))
allocate(wfmt2(ntp,nmtloc,nspinor))
allocate(wfir2(ngrloc,nspinor))
allocate(exmt_(ntp,nmtloc))
allocate(exir_(ngrloc))
allocate(ecmt_(ntp,nmtloc))
allocate(ecir_(ngrloc))
allocate(vxmt_(ntp,nmtloc,nspinor))
allocate(vxir_(ngrloc,nspinor))
allocate(vcmt_(ntp,nmtloc,nspinor))
allocate(vcir_(ngrloc,nspinor))
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  if (sic_apply(n).eq.2) then
    excwanmt=0.d0
    excwanir=0.d0
    do it=1,sic_orbitals%ntr
!-----------------!
! muffin-tin part !
!-----------------!
      if (sic_orbitals%twanmtuc(it,n)) then
        wfmt=zzero
        vxmt_=0.d0
        vcmt_=0.d0
        exmt_=0.d0
        ecmt_=0.d0
        do ispn=1,nspinor
! wfmt(tp,r)=\sum_{lm} f_{lm}(r) * Y_{lm}(tp)
          call zgemm('T','N',ntp,nmtloc,lmmaxvr,zone,ylm,lmmaxvr,&
            sic_orbitals%wanmt(1,1,it,ispn,j),lmmaxvr,zzero,wfmt,ntp)
! rho(tp,r)=|wf(tp,r)|^2
          wfmt2(:,:,ispn)=dreal(dconjg(wfmt(:,:))*wfmt(:,:))
        enddo
! compute XC potential and energy density
        if (spinpol) then
          call xcifc(xctype,n=ntp*nmtloc,rhoup=wfmt2(1,1,1),rhodn=wfmt2(1,1,2),&
            ex=exmt_,ec=ecmt_,vxup=vxmt_(1,1,1),vxdn=vxmt_(1,1,2),&
            vcup=vcmt_(1,1,1),vcdn=vcmt_(1,1,2))
        else
          call xcifc(xctype,n=ntp*nmtloc,rho=wfmt2,ex=exmt_,ec=ecmt_,vx=vxmt_,&
            vc=vcmt_)
        endif
! save total XC potential to vxmt_ 
        do ispn=1,nspinor
          vxmt_(:,:,ispn)=vxmt_(:,:,ispn)+vcmt_(:,:,ispn)
        enddo
! expand XC potential in real spherical harmonics and add it to vhxcmt
        do ispn=1,nspinor
          call dgemm('T','N',lmmaxvr,nmtloc,ntp,1.d0,rlmb,ntp,&
            vxmt_(1,1,ispn),ntp,1.d0,vhxcmt(1,1,it,ispn,j),lmmaxvr)
        enddo
! save total XC energy to exmt_
        exmt_(:,:)=exmt_(:,:)+ecmt_(:,:)
! expand XC energy in real spherical harmonics
        call dgemm('T','N',lmmaxvr,nmtloc,ntp,1.d0,rlmb,ntp,exmt_,ntp,0.d0,&
          excwanmt(1,1,it),lmmaxvr)
      endif
!-------------------!
! interstitial part !
!-------------------!
      do ispn=1,nspinor
        wfir2(:,ispn)=dreal(dconjg(sic_orbitals%wanir(:,it,ispn,j))*&
          sic_orbitals%wanir(:,it,ispn,j))
      enddo
      ecir_=0.d0
      exir_=0.d0
      vxir_=0.d0
      vcir_=0.d0
      if (spinpol) then
        call xcifc(xctype,n=ngrloc,rhoup=wfir2(:,1),rhodn=wfir2(:,2),ex=exir_,&
          ec=ecir_,vxup=vxir_(:,1),vxdn=vxir_(:,2),vcup=vcir_(:,1),vcdn=vcir_(:,2))
      else
        call xcifc(xctype,n=ngrloc,rho=wfir2,ex=exir_,ec=ecir_,vx=vxir_,vc=vcir_)
      endif
      do ispn=1,nspinor
        vhxcir(:,it,ispn,j)=vhxcir(:,it,ispn,j)+vxir_(:,ispn)+vcir_(:,ispn)
      enddo
      excwanir(:,it)=exir_(:)+ecir_(:)
    enddo !it
    do ispn=1,nspinor
      ene(4,j)=ene(4,j)+sic_int_zdz(sic_orbitals%wanmt(1,1,1,ispn,j),&
        sic_orbitals%wanir(1,1,ispn,j),excwanmt,excwanir,&
        sic_orbitals%wanmt(1,1,1,ispn,j),sic_orbitals%wanir(1,1,ispn,j),&
        sic_orbitals%twanmtuc(1,n))
    enddo
  endif
enddo !n
deallocate(ylm)
deallocate(rlmb)
deallocate(rlmz)
deallocate(wfmt)
deallocate(wfmt2)
deallocate(wfir2)
deallocate(exmt_)
deallocate(exir_)
deallocate(ecmt_)
deallocate(ecir_)
deallocate(vxmt_)
deallocate(vxir_)
deallocate(vcmt_)
deallocate(vcir_)
deallocate(excwanmt)
deallocate(excwanir)
return
end
