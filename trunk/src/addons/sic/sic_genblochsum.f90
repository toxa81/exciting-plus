subroutine sic_genblochsum(twk,twvk)
use modmain
use mod_sic
use mod_ws
implicit none
!
logical, intent(in) :: twk
logical, intent(in) :: twvk
!
integer ik,ikloc,j,n,it,ias,is,ia,ir,itp,ispn,ntrloc,itloc,ig
integer irloc,nrloc
real(8) x(3),x0(3)
complex(8) zt1(nspinor),zt2(nspinor),z1,expikt
real(8), allocatable :: vtcmt(:,:,:,:)
real(8), allocatable :: vtcir(:,:)
complex(8), allocatable :: f1tp(:,:,:,:)
complex(8), allocatable :: f1ktp(:,:,:,:,:)
complex(8), allocatable :: f2tp(:,:,:,:)
complex(8), allocatable :: f2ktp(:,:,:,:,:)
complex(8), allocatable :: f1ir(:,:)
complex(8), allocatable :: f1kir(:,:,:)
complex(8), allocatable :: f2ir(:,:)
complex(8), allocatable :: f2kir(:,:,:)
real(8) t0,t1
!
if (twk) then
  s_wkmt=zzero
  s_wkit=zzero
endif
if (twvk) then
  s_wvkmt=zzero
  s_wvkit=zzero
endif
if (.not.tsic_wv) return
call timer_start(90,reset=.true.)
call timer_reset(91)
call timer_reset(92)
call timer_reset(93)
ntrloc=mpi_grid_map(sic_ntr,dim2)
!n=mpi_grid_map(nkpt,dim_k,offs=koffs)
! allocate arrays
allocate(vtcmt(3,mt_ntp,nrmtmax,natmtot))
allocate(vtcir(3,ngrtot))
if (twk) then
  allocate(f1tp(mt_ntp,nrmtmax,natmtot,nspinor))
  allocate(f1ktp(mt_ntp,nrmtmax,natmtot,nspinor,nkptloc))
  allocate(f1ir(ngrtot,nspinor))
  allocate(f1kir(ngrtot,nspinor,nkptloc))
endif
if (twvk) then
  allocate(f2tp(mt_ntp,nrmtmax,natmtot,nspinor))
  allocate(f2ktp(mt_ntp,nrmtmax,natmtot,nspinor,nkptloc))
  allocate(f2ir(ngrtot,nspinor))
  allocate(f2kir(ngrtot,nspinor,nkptloc))
endif
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  if (twk) then
    f1ktp=zzero
    f1kir=zzero
  endif
  if (twvk) then
    f2ktp=zzero
    f2kir=zzero
  endif
! loop over translations: fk(r)=\sum_{T}e^{ikT}f(r-T)=\sum_{T}e^{-ikT}f(r+T) 
  do itloc=1,ntrloc  
    it=mpi_grid_map(sic_ntr,dim2,loc=itloc)
! muffin-tin part
    call timer_start(92)
    vtcmt=0.d0
    if (twk) f1tp=zzero
    if (twvk) f2tp=zzero
    do ias=1,natmtot
      is=ias2is(ias)
      ia=ias2ia(ias)
      nrloc=mpi_grid_map(nrmt(is),dim_k)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(irloc,ir,itp,x0,x,zt1,zt2)
      do irloc=1,nrloc
        ir=mpi_grid_map(nrmt(is),dim_k,loc=irloc)
        do itp=1,mt_ntp
          x0(:)=mt_spx(:,itp)*spr(ir,is)+atposc(:,ia,is)-wanpos(:,n)
          x(:)=x0(:)+sic_vtc(:,it)
          call ws_reduce(x,sic_wan_rwsmax)
          vtcmt(:,itp,ir,ias)=x(:)-x0(:)  !TODO: check + or - ; check t= argument
          if (twk.and.twvk) then
            call s_spinor_func_val(x,s_wlm(1,1,1,j),zt1,&
              &s_wvlm(1,1,1,j),zt2,rcutoff=s_rmax)
            f1tp(itp,ir,ias,:)=zt1(:)
            f2tp(itp,ir,ias,:)=zt2(:)
          else
            if (twk) then
              call s_spinor_func_val(x,s_wlm(1,1,1,j),zt1,rcutoff=s_rmax)
              f1tp(itp,ir,ias,:)=zt1(:)
            endif
            if (twvk) then
              call s_spinor_func_val(x,s_wvlm(1,1,1,j),zt2,rcutoff=s_rmax)
              f2tp(itp,ir,ias,:)=zt2(:)
            endif
          endif
        enddo !itp
      enddo !irloc
!$OMP END PARALLEL DO
    enddo !ias
    call mpi_grid_reduce(vtcmt(1,1,1,1),3*mt_ntp*nrmtmax*natmtot,&
      &dims=(/dim_k/),all=.true.)
    if (twk) call mpi_grid_reduce(f1tp(1,1,1,1),mt_ntp*nrmtmax*natmtot*nspinor,&
      &dims=(/dim_k/),all=.true.)
    if (twvk) call mpi_grid_reduce(f2tp(1,1,1,1),mt_ntp*nrmtmax*natmtot*nspinor,&
      &dims=(/dim_k/),all=.true.)
    call timer_stop(92)
! multiply by e^{-ikT} and add to Bloch-sum 
    call timer_start(93)
    t0=0
    expikt=zone
    do ikloc=1,nkptloc
      ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
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
              if (twk) f1ktp(itp,ir,ias,ispn,ikloc)=&
                &f1ktp(itp,ir,ias,ispn,ikloc)+f1tp(itp,ir,ias,ispn)*expikt
              if (twvk) f2ktp(itp,ir,ias,ispn,ikloc)=&
                &f2ktp(itp,ir,ias,ispn,ikloc)+f2tp(itp,ir,ias,ispn)*expikt
            enddo !ispn
          enddo !itp
        enddo !ir
      enddo !ias
    enddo !ikloc
    call timer_stop(93)
! interstitial part
    call timer_start(92)
    vtcir=0.d0
    if (twk) f1ir=zzero
    if (twvk) f2ir=zzero
    nrloc=mpi_grid_map(ngrtot,dim_k)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(irloc,ir,x0,x,zt1,zt2)
    do irloc=1,nrloc
      ir=mpi_grid_map(ngrtot,dim_k,loc=irloc)
      x0(:)=vgrc(:,ir)-wanpos(:,n)
      x(:)=x0(:)+sic_vtc(:,it)
      call ws_reduce(x,sic_wan_rwsmax)
      vtcir(:,ir)=x(:)-x0(:)
      if (twk.and.twvk) then
        call s_spinor_func_val(x,s_wlm(1,1,1,j),zt1,&
          &s_wvlm(1,1,1,j),zt2,rcutoff=s_rmax)
        f1ir(ir,:)=zt1(:)
        f2ir(ir,:)=zt2(:)
      else
        if (twk) then
          call s_spinor_func_val(x,s_wlm(1,1,1,j),zt1,rcutoff=s_rmax)
          f1ir(ir,:)=zt1(:)
        endif
        if (twvk) then
          call s_spinor_func_val(x,s_wvlm(1,1,1,j),zt2,rcutoff=s_rmax)
          f2ir(ir,:)=zt2(:)
        endif
      endif
    enddo !irloc
!$OMP END PARALLEL DO
    call mpi_grid_reduce(vtcir(1,1),3*ngrtot,dims=(/dim_k/),all=.true.)
    if (twk) call mpi_grid_reduce(f1ir(1,1),ngrtot*nspinor,&
      &dims=(/dim_k/),all=.true.)
    if (twvk) call mpi_grid_reduce(f2ir(1,1),ngrtot*nspinor,&
      &dims=(/dim_k/),all=.true.)
    call timer_stop(92)
! multiply by e^{-ikT} and add to Bloch-sum 
    call timer_start(93)
    t0=0
    expikt=zone
    do ikloc=1,nkptloc
      ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
      do ir=1,ngrtot
        t1=dot_product(vkc(:,ik),vtcir(:,ir))
        if (abs(t1-t0).gt.1d-10) then
          expikt=exp(-zi*t1)
          t0=t1
        endif
        do ispn=1,nspinor
          if (twk) f1kir(ir,ispn,ikloc)=f1kir(ir,ispn,ikloc)+f1ir(ir,ispn)*expikt
          if (twvk) f2kir(ir,ispn,ikloc)=f2kir(ir,ispn,ikloc)+f2ir(ir,ispn)*expikt
        enddo !ispn
      enddo !ikloc
    enddo !ir
  enddo !itloc
  call timer_stop(93)
  do ikloc=1,nkptloc
    ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
    if (twk) then
      call mpi_grid_reduce(f1ktp(1,1,1,1,ikloc),&
        &mt_ntp*nrmtmax*natmtot*nspinor,dims=(/dim2/),all=.true.)
      call mpi_grid_reduce(f1kir(1,1,ikloc),ngrtot*nspinor,&
        &dims=(/dim2/),all=.true.)
    endif
    if (twvk) then
      call mpi_grid_reduce(f2ktp(1,1,1,1,ikloc),&
        &mt_ntp*nrmtmax*natmtot*nspinor,dims=(/dim2/),all=.true.)
      call mpi_grid_reduce(f2kir(1,1,ikloc),ngrtot*nspinor,&
        &dims=(/dim2/),all=.true.)    
    endif
    do ispn=1,nspinor
! convert to lm expansion
      do ias=1,natmtot
        if (twk) call zgemm('T','N',nrmtmax,lmmaxapw,mt_ntp,zone,&
          &f1ktp(1,1,ias,ispn,ikloc),mt_ntp,mt_ylmb,mt_ntp,&
          &zzero,s_wkmt(1,1,ias,ispn,j,ikloc),nrmtmax)
        if (twvk) call zgemm('T','N',nrmtmax,lmmaxapw,mt_ntp,zone,&
          &f2ktp(1,1,ias,ispn,ikloc),mt_ntp,mt_ylmb,mt_ntp,&
          &zzero,s_wvkmt(1,1,ias,ispn,j,ikloc),nrmtmax)
      enddo
! convert to plane-wave expansion
! we want to compute <G+k|w_{nk}> where |G+k> = Omega^{-1/2}*exp^{i(G+k)r}
! straightforward way:
!   do ir=1,ngrtot
!    expikr=exp(-zi*dot_product(vgkc(:,ig,1,ikloc),vgrc(:,ir)))/sqrt(omega)
!    z1=z1+expikr*wkir(ir,ispn,ikloc)*cfunir(ir)*omega/dble(ngrtot)
!   enddo
! but this can be done much faster with FFT
      do ir=1,ngrtot
! remove exp^{-ikr} phase, multiply by step function and by sqrt(omega)
        z1=exp(-zi*dot_product(vkc(:,ik),vgrc(:,ir)))*sqrt(omega)
        if (twk) f1kir(ir,ispn,ikloc)=f1kir(ir,ispn,ikloc)*cfunir(ir)*z1 !TODO: check if cfunir is really necessary
        if (twvk) f2kir(ir,ispn,ikloc)=f2kir(ir,ispn,ikloc)*cfunir(ir)*z1
      enddo !ir
      if (twk) call zfftifc(3,ngrid,-1,f1kir(1,ispn,ikloc))
      if (twvk) call zfftifc(3,ngrid,-1,f2kir(1,ispn,ikloc))
      do ig=1,ngk(1,ik)
        if (twk) s_wkit(ig,ispn,j,ikloc)=f1kir(igfft(igkig(ig,1,ikloc)),ispn,ikloc)
        if (twvk) s_wvkit(ig,ispn,j,ikloc)=f2kir(igfft(igkig(ig,1,ikloc)),ispn,ikloc)
      enddo !ig
    enddo !ispn
  enddo !ikloc
enddo !j
deallocate(vtcmt)
deallocate(vtcir)
if (twk) then
  deallocate(f1tp)
  deallocate(f1ktp)
  deallocate(f1ir)
  deallocate(f1kir)
endif
if (twvk) then
  deallocate(f2tp)
  deallocate(f2ktp)
  deallocate(f2ir)
  deallocate(f2kir)
endif
call timer_stop(90)
if (mpi_grid_root()) then
  write(*,'("[sic_genblochsum] total time   : ",F12.4," sec.")')timer_get_value(90)
  write(*,'("[sic_genblochsum] wf(r)        : ",F12.4," sec.")')timer_get_value(92)
  write(*,'("[sic_genblochsum] summation    : ",F12.4," sec.")')timer_get_value(93)
endif
return
end
