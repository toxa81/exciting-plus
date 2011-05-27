subroutine sic_genblochsum
use modmain
use mod_sic
use mod_ws
implicit none
integer ik,ikloc,j,n,it,ias,is,ia,ir,itp,ispn,ntrloc,itloc,ig,koffs
integer irloc,nrloc
real(8) x(3),vtc(3),x0(3)
complex(8) zt1(nspinor),zt2(nspinor),z1,z2,expikt,expikr
complex(8), allocatable :: wkir(:,:,:)
complex(8), allocatable :: wvkir(:,:,:)
real(8), allocatable :: vtcmt(:,:,:,:)
real(8), allocatable :: vtcir(:,:)
complex(8), allocatable :: f1mt(:,:,:,:)
complex(8), allocatable :: f2mt(:,:,:,:)
complex(8), allocatable :: f1ir(:,:)
complex(8), allocatable :: f2ir(:,:)
real(8) t0,t1
!
s_wkmt=zzero
s_wvkmt=zzero
s_wkit=zzero
s_wvkit=zzero
if (.not.tsic_wv) return
call timer_start(90,reset=.true.)
call timer_reset(91)
call timer_reset(92)
call timer_reset(93)
ntrloc=mpi_grid_map(sic_blochsum%ntr,dim2)
n=mpi_grid_map(nkpt,dim_k,offs=koffs)
allocate(wkir(ngrtot,nspinor,nkptloc))
allocate(wvkir(ngrtot,nspinor,nkptloc))

allocate(vtcmt(3,mt_ntp,nrmtmax,natmtot))
allocate(vtcir(3,ngrtot))
allocate(f1mt(mt_ntp,nrmtmax,natmtot,nspinor))
allocate(f2mt(mt_ntp,nrmtmax,natmtot,nspinor))
allocate(f1ir(ngrtot,nspinor))
allocate(f2ir(ngrtot,nspinor))

do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  wkir=zzero
  wvkir=zzero
  do itloc=1,ntrloc
    call timer_start(92)
    it=mpi_grid_map(sic_blochsum%ntr,dim2,loc=itloc)
    vtcmt=0.d0
    vtcir=0.d0
    f1mt=zzero
    f2mt=zzero
    f1ir=zzero
    f2ir=zzero
! muffin-tin part
    do ias=1,natmtot
      is=ias2is(ias)
      ia=ias2ia(ias)
      nrloc=mpi_grid_map(nrmt(is),dim_k)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(irloc,ir,itp,x0,x,zt1,zt2)
      do irloc=1,nrloc
        ir=mpi_grid_map(nrmt(is),dim_k,loc=irloc)
        do itp=1,mt_ntp
          x0(:)=mt_spx(:,itp)*spr(ir,is)+atposc(:,ia,is)-wanpos(:,n)
          x(:)=x0(:)+sic_blochsum%vtc(:,it)
          call ws_reduce(x,sic_wan_rwsmax)
          vtcmt(:,itp,ir,ias)=x(:)-x0(:)
          call s_spinor_func_val(x,s_wlm(1,1,1,j),zt1,&
            s_wvlm(1,1,1,j),zt2,rcutoff=s_rmax)
          f1mt(itp,ir,ias,:)=zt1(:)
          f2mt(itp,ir,ias,:)=zt2(:)
        enddo
      enddo
!$OMP END PARALLEL DO
    enddo
    call mpi_grid_reduce(vtcmt(1,1,1,1),3*mt_ntp*nrmtmax*natmtot,&
      dims=(/dim_k/),all=.true.)
    call mpi_grid_reduce(f1mt(1,1,1,1),mt_ntp*nrmtmax*natmtot*nspinor,&
      dims=(/dim_k/),all=.true.)
    call mpi_grid_reduce(f2mt(1,1,1,1),mt_ntp*nrmtmax*natmtot*nspinor,&
      dims=(/dim_k/),all=.true.)
    call timer_stop(92)
    call timer_start(93)
    t0=0
    expikt=zone
    do ikloc=1,nkptloc
      ik=ikloc+koffs
      do ias=1,natmtot
        is=ias2is(ias)
        ia=ias2ia(ias)
        do ir=1,nrmt(is)
          do itp=1,mt_ntp
            t1=dot_product(vkc(:,ik),vtcmt(:,itp,ir,ias))
            if (abs(t1-t0).gt.1d-10) then
              expikt=exp(-zi*t1)
              t0=t1
            endif
            do ispn=1,nspinor
              s_wkmt(itp,ir,ias,ispn,j,ikloc)=&
                s_wkmt(itp,ir,ias,ispn,j,ikloc)+f1mt(itp,ir,ias,ispn)*expikt
              s_wvkmt(itp,ir,ias,ispn,j,ikloc)=&
                s_wvkmt(itp,ir,ias,ispn,j,ikloc)+f2mt(itp,ir,ias,ispn)*expikt
            enddo !ispn
          enddo !
        enddo !
      enddo !
    enddo !
    call timer_stop(93)
    call timer_start(92)
! interstitial part
    nrloc=mpi_grid_map(ngrtot,dim_k)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(irloc,ir,x0,x,zt1,zt2)
    do irloc=1,nrloc
      ir=mpi_grid_map(ngrtot,dim_k,loc=irloc)
      x0(:)=vgrc(:,ir)-wanpos(:,n)
      x(:)=x0(:)+sic_blochsum%vtc(:,it)
      call ws_reduce(x,sic_wan_rwsmax)
      vtcir(:,ir)=x(:)-x0(:)
      call s_spinor_func_val(x,s_wlm(1,1,1,j),zt1,&
        s_wvlm(1,1,1,j),zt2,rcutoff=s_rmax)
      f1ir(ir,:)=zt1(:)
      f2ir(ir,:)=zt2(:)
    enddo
!$OMP END PARALLEL DO
    call mpi_grid_reduce(vtcir(1,1),3*ngrtot,dims=(/dim_k/),all=.true.)
    call mpi_grid_reduce(f1ir(1,1),ngrtot*nspinor,dims=(/dim_k/),all=.true.)
    call mpi_grid_reduce(f2ir(1,1),ngrtot*nspinor,dims=(/dim_k/),all=.true.)
    call timer_stop(92)
    call timer_start(93)
    t0=0
    expikt=zone
    do ikloc=1,nkptloc
      ik=ikloc+koffs
      do ir=1,ngrtot
        t1=dot_product(vkc(:,ik),vtcir(:,ir))
        if (abs(t1-t0).gt.1d-10) then
          expikt=exp(-zi*t1)
          t0=t1
        endif
        do ispn=1,nspinor
          wkir(ir,ispn,ikloc)=wkir(ir,ispn,ikloc)+f1ir(ir,ispn)*expikt
          wvkir(ir,ispn,ikloc)=wvkir(ir,ispn,ikloc)+f2ir(ir,ispn)*expikt
        enddo !ispn
      enddo !ikloc
    enddo !ir
  enddo !itloc
  call timer_stop(93)
  do ikloc=1,nkptloc
    call mpi_grid_reduce(s_wkmt(1,1,1,1,j,ikloc),&
      mt_ntp*nrmtmax*natmtot*nspinor,dims=(/dim2/),all=.true.)
    call mpi_grid_reduce(s_wvkmt(1,1,1,1,j,ikloc),&
      mt_ntp*nrmtmax*natmtot*nspinor,dims=(/dim2/),all=.true.)
    call mpi_grid_reduce(wkir(1,1,ikloc),ngrtot*nspinor,&
      dims=(/dim2/),all=.true.)
    call mpi_grid_reduce(wvkir(1,1,ikloc),ngrtot*nspinor,&
      dims=(/dim2/),all=.true.)
  enddo
! convert to plane-wave expansion
! we want to compute <G+k|w_{nk}> where |G+k> = Omega^{-1/2}*exp^{i(G+k)r}
! straightforward way:
! do ir=1,ngrtot
!  expikr=exp(-zi*dot_product(vgkc(:,ig,1,ikloc),vgrc(:,ir)))/sqrt(omega)
!  z1=z1+expikr*wkir(ir,ispn,ikloc)*cfunir(ir)*omega/dble(ngrtot)
! enddo
! but this can be done much faster with FFT
  call timer_start(91)
  do ispn=1,nspinor
    do ikloc=1,nkptloc
      ik=ikloc+koffs
! remove exp^{-ikr} phase, multiply by step function and by sqrt(omega)
      do ir=1,ngrtot
        z1=exp(-zi*dot_product(vkc(:,ik),vgrc(:,ir)))*sqrt(omega)
        wkir(ir,ispn,ikloc)=wkir(ir,ispn,ikloc)*cfunir(ir)*z1
        wvkir(ir,ispn,ikloc)=wvkir(ir,ispn,ikloc)*cfunir(ir)*z1
      enddo
      call zfftifc(3,ngrid,-1,wkir(1,ispn,ikloc))
      call zfftifc(3,ngrid,-1,wvkir(1,ispn,ikloc))
      do ig=1,ngk(1,ik)
        s_wkit(ig,ispn,j,ikloc)=wkir(igfft(igkig(ig,1,ikloc)),ispn,ikloc)
        s_wvkit(ig,ispn,j,ikloc)=wvkir(igfft(igkig(ig,1,ikloc)),ispn,ikloc)
      enddo
      !do ig=1,ngk(1,ik)
      !  z1=zzero
      !  z2=zzero
      !  do ir=1,ngrtot
      !    expikr=exp(-zi*dot_product(vgkc(:,ig,1,ikloc),vgrc(:,ir)))/sqrt(omega)
      !    z1=z1+expikr*wkir(ir,ispn,ikloc)*cfunir(ir)*omega/dble(ngrtot)
      !    z2=z2+expikr*wvkir(ir,ispn,ikloc)*cfunir(ir)*omega/dble(ngrtot)
      !  enddo
      !  s_wkit(ig,ispn,j,ikloc)=z1
      !  s_wvkit(ig,ispn,j,ikloc)=z2
      !enddo !ig
    enddo !ikloc
  enddo !ispn
  call timer_stop(91)
enddo !j
deallocate(wkir,wvkir)
deallocate(vtcmt,vtcir,f1mt,f2mt,f1ir,f2ir)
call timer_stop(90)
if (mpi_grid_root()) then
  write(*,'("[sic_blochsum] total time   : ",F12.4," sec.")')timer_get_value(90)
  write(*,'("[sic_blochsum] wf(r)        : ",F12.4," sec.")')timer_get_value(92)
  write(*,'("[sic_blochsum] summation    : ",F12.4," sec.")')timer_get_value(93)
  write(*,'("[sic_blochsum] pw expansion : ",F12.4," sec.")')timer_get_value(91)
endif
return
end
