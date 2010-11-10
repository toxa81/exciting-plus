subroutine sic_genvxc(vhxcmt,vhxcir,ene)
use modmain
use mod_lf
use modxcifc
implicit none
real(8), intent(out) :: vhxcmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor,nwannloc)
real(8), intent(out) :: vhxcir(ngrtot,ntrloc,nspinor,nwannloc)
complex(8), intent(out) :: ene(4,nwann)
integer ntp,itp,lm,n,itloc,ias,ispn,nloc,it,l
real(8), allocatable :: tp(:,:)
complex(8), allocatable :: ylm(:,:)
complex(8), allocatable :: rlmz(:,:)
real(8) rlm(lmmaxvr),theta,phi
real(8), allocatable :: rlmb(:,:)
real(8), allocatable :: vxcwanmt(:,:,:,:,:)
real(8), allocatable :: vxcwanir(:,:,:)
real(8), allocatable :: excwanmt(:,:,:,:)
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
real(8), allocatable :: vxcmt_(:,:,:)
real(8), allocatable :: excmt_(:,:)
real(8), allocatable :: x(:),w(:),wtp(:)
integer nt,np,ip,i,j
real(8), allocatable :: spx(:,:)
real(8) t1

! XC part of Wannier function potential
allocate(vxcwanmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor))
allocate(vxcwanir(ngrtot,ntrloc,nspinor))
! XC energy density of Wannier function
allocate(excwanmt(lmmaxvr,nrmtmax,natmtot,ntrloc))
allocate(excwanir(ngrtot,ntrloc))



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

allocate(wfmt(ntp,nrmtmax))
allocate(wfmt2(ntp,nrmtmax,nspinor))
allocate(wfir2(ngrtot,nspinor))
allocate(exmt_(ntp,nrmtmax))
allocate(exir_(ngrtot))
allocate(ecmt_(ntp,nrmtmax))
allocate(ecir_(ngrtot))
allocate(excmt_(ntp,nrmtmax))
allocate(vxmt_(ntp,nrmtmax,nspinor))
allocate(vxir_(ngrtot,nspinor))
allocate(vcmt_(ntp,nrmtmax,nspinor))
allocate(vcir_(ngrtot,nspinor))
allocate(vxcmt_(ntp,nrmtmax,nspinor))
do nloc=1,nwannloc
  n=mpi_grid_map(nwann,dim_k,loc=nloc)
  vxcwanmt=0.d0
  vxcwanir=0.d0
  excwanmt=0.d0
  excwanir=0.d0
  do itloc=1,ntrloc
    it=mpi_grid_map(ntr,dim_t,loc=itloc)
!-----------------!
! muffin-tin part !
!-----------------!
    do ias=1,natmtot
      if (twanmt(ias,it,n)) then
        wfmt=zzero
        vxmt_=0.d0
        vcmt_=0.d0
        exmt_=0.d0
        ecmt_=0.d0
! compute charge density on a sphere
!   rho(tp,r)=|wf(tp,r)|^2
!   wf(tp,r)=\sum_{lm} f_{lm}(r) * Y_{lm}(tp)
        do ispn=1,nspinor
          call zgemm('T','N',ntp,nrmt(ias2is(ias)),lmmaxvr,zone,ylm,lmmaxvr,&
            wanmt(1,1,ias,itloc,ispn,nloc),lmmaxvr,zzero,wfmt,ntp)
          wfmt2(:,:,ispn)=dreal(dconjg(wfmt(:,:))*wfmt(:,:))
        enddo
! compute XC potential and energy density
        if (spinpol) then
          call xcifc(xctype,n=ntp*nrmtmax,rhoup=wfmt2(1,1,1),rhodn=wfmt2(1,1,2),&
            ex=exmt_,ec=ecmt_,vxup=vxmt_(1,1,1),vxdn=vxmt_(1,1,2),vcup=vcmt_(1,1,1),&
            vcdn=vcmt_(1,1,2))
       else
          call xcifc(xctype,n=ntp*nrmtmax,rho=wfmt2,ex=exmt_,ec=ecmt_,vx=vxmt_,vc=vcmt_)
        endif
! save XC potential
        do ispn=1,nspinor
          vxcmt_(:,:,ispn)=vxmt_(:,:,ispn)+vcmt_(:,:,ispn)
        enddo
! save XC energy
        excmt_(:,:)=exmt_(:,:)+ecmt_(:,:)
! expand XC potential in spherical harmonics
        do ispn=1,nspinor
          call dgemm('T','N',lmmaxvr,nrmt(ias2is(ias)),ntp,1.d0,rlmb,ntp,&
            vxcmt_(1,1,ispn),ntp,0.d0,vxcwanmt(1,1,ias,itloc,ispn),lmmaxvr)
        enddo
! expand XC energy in spherical harmonics
        call dgemm('T','N',lmmaxvr,nrmt(ias2is(ias)),ntp,1.d0,rlmb,ntp,&
          excmt_,ntp,0.d0,excwanmt(1,1,ias,itloc),lmmaxvr)
      endif
    enddo !ias
!-------------------!
! interstitial part !
!-------------------!
    do ispn=1,nspinor
      wfir2(:,ispn)=dreal(dconjg(wanir(:,itloc,ispn,nloc))*wanir(:,itloc,ispn,nloc))
    enddo
    ecir_=0.d0
    exir_=0.d0
    vxir_=0.d0
    vcir_=0.d0
    if (spinpol) then
      call xcifc(xctype,n=ngrtot,rhoup=wfir2(:,1),rhodn=wfir2(:,2),ex=exir_,&
        ec=ecir_,vxup=vxir_(:,1),vxdn=vxir_(:,2),vcup=vcir_(:,1),vcdn=vcir_(:,2))
    else
      call xcifc(xctype,n=ngrtot,rho=wfir2,ex=exir_,ec=ecir_,vx=vxir_,vc=vcir_)
    endif
    do ispn=1,nspinor
      vxcwanir(:,itloc,ispn)=vxir_(:,ispn)+vcir_(:,ispn)
    enddo
    excwanir(:,itloc)=exir_(:)+ecir_(:)
  enddo !itloc
  do ispn=1,nspinor
    vhxcmt(:,:,:,:,ispn,nloc)=vhxcmt(:,:,:,:,ispn,nloc)+vxcwanmt(:,:,:,:,ispn)
    vhxcir(:,:,ispn,nloc)=vhxcir(:,:,ispn,nloc)+vxcwanir(:,:,ispn)
  enddo
  do ispn=1,nspinor
    ene(2,n)=ene(2,n)+lf_intgr_zdz(wanmt(1,1,1,1,ispn,nloc),wanir(1,1,ispn,nloc),&
      excwanmt,excwanir,(/0,0,0/),wanmt(1,1,1,1,ispn,nloc),wanir(1,1,ispn,nloc))
  enddo
enddo !nloc
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
deallocate(excmt_)
deallocate(vxmt_)
deallocate(vxir_)
deallocate(vcmt_)
deallocate(vcir_)
deallocate(vxcmt_)
deallocate(vxcwanmt)
deallocate(vxcwanir)
deallocate(excwanmt)
deallocate(excwanir)
return
end