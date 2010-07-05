subroutine sic_genvxc(exc)
use modmain
use mod_lf
use modxcifc
implicit none
complex(8), intent(out) :: exc(nwann)
integer ntp,itp,lm,n,m,itloc,irloc,nrmtloc,ias,ispn,ir
complex(8) zt1
real(8), allocatable :: tp(:,:)
real(8), allocatable :: vx0(:,:),vc0(:,:),ex0(:),ec0(:),rho0(:,:)
complex(8), allocatable :: zrho0(:,:),zvxc(:,:),zexc(:)
complex(8), allocatable :: ylm(:,:)
complex(8), allocatable :: ylmc(:,:)
complex(8), allocatable :: vxcwanmt(:,:,:,:,:)
complex(8), allocatable :: vxcwanir(:,:,:)
complex(8), allocatable :: excwanmt(:,:,:,:,:)
complex(8), allocatable :: excwanir(:,:,:)
real(8), external :: gaunt

! XC part of Wannier function potential
allocate(vxcwanmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor))
allocate(vxcwanir(ngrtot,ntrloc,nspinor))
! XC energy density of Wannier function
allocate(excwanmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor))
allocate(excwanir(ngrtot,ntrloc,nspinor))

! make dens mesh of (theta,phi) points on the sphere
ntp=1000
allocate(tp(2,ntp))
allocate(ylm(lmmaxvr,ntp))
allocate(ylmc(ntp,lmmaxvr))
call sphcover(ntp,tp)
do itp=1,ntp 
  call genylm(lmaxvr,tp(1,itp),ylm(1,itp))
  do lm=1,lmmaxvr
    ylmc(itp,lm)=dconjg(ylm(lm,itp))
  enddo
enddo
m=max(ntp,ngrtot)
allocate(rho0(m,nspinor),vx0(m,nspinor),vc0(m,nspinor))
allocate(ex0(m),ec0(m))
allocate(zrho0(m,nspinor))
allocate(zvxc(m,nspinor))
allocate(zexc(m))
do n=1,nwann
  vxcwanmt=zzero
  vxcwanir=zzero
  excwanmt=zzero
  excwanir=zzero
  do itloc=1,ntrloc
! muffin-tin part
    do ias=1,natmtot
      nrmtloc=mpi_grid_map(nrmt(ias2is(ias)),dim_k)
      do irloc=1,nrmtloc
        ir=mpi_grid_map(nrmt(ias2is(ias)),dim_k,loc=irloc)
! compute charge density on a sphere        
        do ispn=1,nspinor
          !do itp=1,ntp
          !  zt1=zzero
          !  do lm=1,lmmaxvr
          !    zt1=zt1+wanmt(lm,ir,ias,itloc,ispn,n)*ylm(lm,itp)
          !  enddo
          !  rho0(itp,ispn)=abs(zt1)**2
          !enddo !itp
          call zgemv('T',lmmaxvr,ntp,zone,ylm,lmmaxvr,&
            wanmt(1,ir,ias,itloc,ispn,n),1,zzero,zrho0(1,ispn),1)
          rho0(1:ntp,ispn)=dreal(dconjg(zrho0(1:ntp,ispn))*zrho0(1:ntp,ispn))
        enddo !ispn
        if (spinpol) then
          call xcifc(xctype,n=ntp,rhoup=rho0(:,1),rhodn=rho0(:,2),ex=ex0(:),&
            ec=ec0(:),vxup=vx0(:,1),vxdn=vx0(:,2),vcup=vc0(:,1),vcdn=vc0(:,2))
        else
          call xcifc(xctype,n=ntp,rho=rho0(:,1),ex=ex0(:),ec=ec0(:),&
            vx=vx0(:,1),vc=vc0(:,1))
        endif
! save XC potential
        do ispn=1,nspinor
          zvxc(1:ntp,ispn)=dcmplx(vx0(1:ntp,ispn)+vc0(1:ntp,ispn),0.d0)
        enddo
! save XC energy
        zexc(1:ntp)=dcmplx(ex0(1:ntp)+ec0(1:ntp),0.d0)
! expand XC potential in spherical harmonics
        do ispn=1,nspinor
          !do lm=1,lmmaxvr
          !  zt1=zzero
          !  do itp=1,ntp
          !    zt1=zt1+dconjg(ylm(lm,itp))*vx0(itp,ispn)
          !  enddo
          !  vxcwanmt(lm,ir,ias,itloc,ispn)=fourpi*zt1/ntp
          !enddo !lm
          call zgemv('T',ntp,lmmaxvr,dcmplx(fourpi/ntp,0.d0),ylmc,ntp,&
            zvxc(1,ispn),1,zzero,vxcwanmt(1,ir,ias,itloc,ispn),1)
        enddo !ispn
! expand XC energy in spherical harmonics
        !do lm=1,lmmaxvr
        !  zt1=zzero
        !  do itp=1,ntp
        !    zt1=zt1+dconjg(ylm(lm,itp))*ex0(itp)
        !  enddo
        !  excwanmt(lm,ir,ias,itloc,1)=fourpi*zt1/ntp
        !enddo !lm
        call zgemv('T',ntp,lmmaxvr,dcmplx(fourpi/ntp,0.d0),ylmc,ntp,&
          zexc,1,zzero,excwanmt(1,ir,ias,itloc,1),1)        
      enddo !irloc
    enddo !ias
    do ispn=1,nspinor
      call mpi_grid_reduce(vxcwanmt(1,1,1,itloc,ispn),lmmaxvr*nrmtmax*natmtot,&
        dims=(/dim_k/),all=.true.)
    enddo
    call mpi_grid_reduce(excwanmt(1,1,1,itloc,1),lmmaxvr*nrmtmax*natmtot,&
      dims=(/dim_k/),all=.true.)
    if (spinpol) excwanmt(:,:,:,itloc,2)=excwanmt(:,:,:,itloc,1)
! interstitial part
    do ispn=1,nspinor
      rho0(1:ngrtot,ispn)=dreal(dconjg(wanir(:,itloc,ispn,n))*wanir(:,itloc,ispn,n))
    enddo
    if (spinpol) then
      call xcifc(xctype,n=ngrtot,rhoup=rho0(:,1),rhodn=rho0(:,2),ex=ex0(:),&
        ec=ec0(:),vxup=vx0(:,1),vxdn=vx0(:,2),vcup=vc0(:,1),vcdn=vc0(:,2))
    else
      call xcifc(xctype,n=ngrtot,rho=rho0(:,1),ex=ex0(:),ec=ec0(:),&
        vx=vx0(:,1),vc=vc0(:,1))
    endif
    do ispn=1,nspinor
      vxcwanir(1:ngrtot,itloc,ispn)=vx0(1:ngrtot,ispn)+vc0(1:ngrtot,ispn)
    enddo
    excwanir(1:ngrtot,itloc,1)=ex0(1:ngrtot)+ec0(1:ngrtot)
    if (spinpol) excwanir(:,itloc,2)=excwanir(:,itloc,1)
  enddo !itloc
  do ispn=1,nspinor
    call lf_prod(zone,vxcwanmt(1,1,1,1,ispn),vxcwanir(1,1,ispn),&
      wanmt(1,1,1,1,ispn,n),wanir(1,1,ispn,n),zone,&
      vwanmt(1,1,1,1,ispn,n),vwanir(1,1,ispn,n))
    call lf_prod(zone,excwanmt(1,1,1,1,ispn),excwanir(1,1,ispn),&
      wanmt(1,1,1,1,ispn,n),wanir(1,1,ispn,n),zzero,&
      excwanmt(1,1,1,1,ispn),excwanir(1,1,ispn))  
  enddo  
  do ispn=1,nspinor
    exc(n)=exc(n)+lf_dotlf(.true.,(/0,0,0/),excwanmt(1,1,1,1,ispn),&
      excwanir(1,1,ispn),wanmt(1,1,1,1,ispn,n),wanir(1,1,ispn,n))
  enddo
enddo
deallocate(rho0,vx0,vc0,vxcwanmt,vxcwanir,excwanmt,excwanir)
deallocate(ex0,ec0,ylm,tp,zrho0,zvxc,zexc,ylmc)
return
end