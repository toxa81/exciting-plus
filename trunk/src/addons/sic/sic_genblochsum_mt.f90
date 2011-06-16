subroutine sic_genblochsum_mt
use modmain
use mod_sic
use mod_ws
implicit none
integer ik,ikloc,j,n,it,ias,is,ia,itp,ispn,ntrloc,itloc,ig,koffs
integer ir,irloc,nrloc
real(8) x(3),vtc(3),x0(3)
complex(8) zt1(nspinor),zt2(nspinor),z1,z2,expikt,expikr
real(8), allocatable :: vtcmt(:,:,:,:)
complex(8), allocatable :: f1mt(:,:,:,:)
complex(8), allocatable :: f2mt(:,:,:,:)
real(8) t0,t1
!
s_wkmt=zzero
s_wvkmt=zzero
if (.not.tsic_wv) return
call timer_start(90,reset=.true.)
call timer_reset(91)
call timer_reset(92)
call timer_reset(93)
ntrloc=mpi_grid_map(sic_blochsum%ntr,dim2)
n=mpi_grid_map(nkpt,dim_k,offs=koffs)

allocate(vtcmt(3,mt_ntp,nrmtmax,natmtot))
allocate(f1mt(mt_ntp,nrmtmax,natmtot,nspinor))
allocate(f2mt(mt_ntp,nrmtmax,natmtot,nspinor))

do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  do itloc=1,ntrloc
    call timer_start(92)
    it=mpi_grid_map(sic_blochsum%ntr,dim2,loc=itloc)
    vtcmt=0.d0
    f1mt=zzero
    f2mt=zzero
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
    enddo ! ikloc
    call timer_stop(93)
  enddo !itloc
  do ikloc=1,nkptloc
    call mpi_grid_reduce(s_wkmt(1,1,1,1,j,ikloc),&
      mt_ntp*nrmtmax*natmtot*nspinor,dims=(/dim2/),all=.true.)
    call mpi_grid_reduce(s_wvkmt(1,1,1,1,j,ikloc),&
      mt_ntp*nrmtmax*natmtot*nspinor,dims=(/dim2/),all=.true.)
  enddo
enddo !j
deallocate(vtcmt,f1mt,f2mt)
call timer_stop(90)
if (mpi_grid_root()) then
  write(*,'("[sic_blochsum] total time   : ",F12.4," sec.")')timer_get_value(90)
  write(*,'("[sic_blochsum] wf(r)        : ",F12.4," sec.")')timer_get_value(92)
  write(*,'("[sic_blochsum] summation    : ",F12.4," sec.")')timer_get_value(93)
endif
return
end
