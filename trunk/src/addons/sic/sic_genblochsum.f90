subroutine sic_genblochsum
use modmain
use mod_sic
use mod_ws
implicit none
integer ik,ikloc,j,n,it,ias,is,ia,ir,itp,ispn,ntrloc,itloc,ig,koffs
real(8) x(3),vtc(3),x0(3)
complex(8) zt1(nspinor),zt2(nspinor),z1,z2,expikt,expikr
complex(8), allocatable :: wkir(:,:,:)
complex(8), allocatable :: wvkir(:,:,:)
!
s_wkmt=zzero
s_wvkmt=zzero
s_wkit=zzero
s_wvkit=zzero
if (.not.tsic_wv) return
call timer_start(90,reset=.true.)
call timer_reset(91)
ntrloc=mpi_grid_map(sic_blochsum%ntr,dim2)
n=mpi_grid_map(nkpt,dim_k,offs=koffs)
allocate(wkir(ngrtot,nspinor,nkptloc))
allocate(wvkir(ngrtot,nspinor,nkptloc))
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  wkir=zzero
  wvkir=zzero
  do itloc=1,ntrloc
    it=mpi_grid_map(sic_blochsum%ntr,dim2,loc=itloc)
! muffin-tin part
    do ias=1,natmtot
      is=ias2is(ias)
      ia=ias2ia(ias)
      do ir=1,nrmt(is)
        do itp=1,mt_ntp
          x0(:)=mt_spx(:,itp)*spr(ir,is)+atposc(:,ia,is)-wanpos(:,n)
          x(:)=x0(:)+sic_blochsum%vtc(:,it)
          call ws_reduce(x,sic_wan_rwsmax)
          vtc(:)=x(:)-x0(:)
          call s_spinor_func_val(x,s_wlm(1,1,1,j),zt1,&
            s_wvlm(1,1,1,j),zt2,rcutoff=s_rmax)
          do ikloc=1,nkptloc
            ik=ikloc+koffs
            expikt=exp(-zi*dot_product(vkc(:,ik),vtc(:))) 
            do ispn=1,nspinor
              s_wkmt(itp,ir,ias,ispn,j,ikloc)=&
                s_wkmt(itp,ir,ias,ispn,j,ikloc)+zt1(ispn)*expikt
              s_wvkmt(itp,ir,ias,ispn,j,ikloc)=&
                s_wvkmt(itp,ir,ias,ispn,j,ikloc)+zt2(ispn)*expikt
            enddo !ispn
          enddo !ikloc
        enddo !itp
      enddo !ir
    enddo !ias
! interstitial part
    do ir=1,ngrtot
      x0(:)=vgrc(:,ir)-wanpos(:,n)
      x(:)=x0(:)+sic_blochsum%vtc(:,it)
      call ws_reduce(x,sic_wan_rwsmax)
      vtc(:)=x(:)-x0(:)
      call s_spinor_func_val(x,s_wlm(1,1,1,j),zt1,&
        s_wvlm(1,1,1,j),zt2,rcutoff=s_rmax)
      do ikloc=1,nkptloc
        ik=ikloc+koffs
        expikt=exp(-zi*dot_product(vkc(:,ik),vtc(:)))
        do ispn=1,nspinor
          wkir(ir,ispn,ikloc)=wkir(ir,ispn,ikloc)+zt1(ispn)*expikt
          wvkir(ir,ispn,ikloc)=wvkir(ir,ispn,ikloc)+zt2(ispn)*expikt
        enddo !ispn
      enddo !ikloc
    enddo !ir
  enddo !itloc
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
  call timer_start(91)
  do ispn=1,nspinor
    do ikloc=1,nkptloc
      ik=ikloc+koffs
      do ig=1,ngk(1,ik)
        z1=zzero
        z2=zzero
        do ir=1,ngrtot
          expikr=exp(-zi*dot_product(vgkc(:,ig,1,ikloc),vgrc(:,ir)))/sqrt(omega)
          z1=z1+expikr*wkir(ir,ispn,ikloc)*cfunir(ir)*omega/dble(ngrtot)
          z2=z2+expikr*wvkir(ir,ispn,ikloc)*cfunir(ir)*omega/dble(ngrtot)
        enddo
        s_wkit(ig,ispn,j,ikloc)=z1
        s_wvkit(ig,ispn,j,ikloc)=z2
      enddo !ig
    enddo !ikloc
  enddo !ispn
  call timer_stop(91)
enddo !j
deallocate(wkir,wvkir)
call timer_stop(90)
if (mpi_grid_root()) then
  write(*,'("[sic_blochsum] total time   : ",F12.4," sec.")')timer_get_value(90)
  write(*,'("[sic_blochsum] pw expansion : ",F12.4," sec.")')timer_get_value(91)
endif
return
end
