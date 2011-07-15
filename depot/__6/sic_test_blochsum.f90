subroutine sic_test_blochsum(itest,fname)
use modmain
use mod_sic
use mod_madness
use mod_ws
implicit none
integer, intent(in) :: itest
character*(*), intent(in) :: fname
integer ik,ikloc,j,n,jas,it,ias,is,ia,ir,itp,ispn,ntrloc,itloc,lm
integer ig,io,ig1,l,itp1
real(8) x(3),t1,tp(2),vg(3)
complex(8) zt1,zt2(2),expikt,wanval(nspinor),wanval1(nspinor)
complex(8), allocatable :: wgk(:,:)
real(8), allocatable :: jl(:,:)
complex(8), allocatable :: zf(:)
complex(8),allocatable :: ylmgk(:)
integer, external :: hash

complex(8), allocatable :: wnklm_exact(:,:,:,:)
complex(8), allocatable :: wnkir_exact(:,:)
complex(8), allocatable :: wnktp_exact(:,:,:,:)
complex(8), allocatable :: wnktp_bt(:,:,:,:)
complex(8), allocatable :: wnkir_bt(:,:)
complex(8), allocatable :: zprodlm_exact(:,:,:)
complex(8), allocatable :: zprodtp_exact(:,:,:)
complex(8), allocatable :: zprodtp_bt(:,:,:)

complex(8), allocatable :: wankmt(:,:,:,:)
complex(8), allocatable :: wankir(:,:)
complex(8), allocatable :: wantp(:,:,:)
complex(8), allocatable :: wanir(:,:)
complex(8), allocatable :: wnkmt(:,:,:,:,:)
complex(8), allocatable :: wnkmt_(:,:)

allocate(wnklm_exact(lmmaxvr,nrmtmax,natmtot,nspinor))
allocate(wnkir_exact(ngrtot,nspinor))
allocate(wnktp_exact(mt_ntp,nrmtmax,natmtot,nspinor))
allocate(wnktp_bt(mt_ntp,nrmtmax,natmtot,nspinor))
allocate(wnkir_bt(ngrtot,nspinor))

allocate(zprodlm_exact(3,sic_wantran%nwan,nkpt))
allocate(zprodtp_exact(3,sic_wantran%nwan,nkpt))
allocate(zprodtp_bt(3,sic_wantran%nwan,nkpt))


zprodlm_exact=zzero
zprodtp_exact=zzero
zprodtp_bt=zzero

ntrloc=mpi_grid_map(sic_blochsum%ntr,dim2)

if (itest.eq.1) then
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    call elk_load_wann_unk(n)
    do ik=1,1 !nkpt
      wnklm_exact=zzero
      do ispn=1,nspinor
        do ias=1,natmtot
          do ir=1,nrmt(ias2is(ias))
            do lm=1,lmmaxvr
              do io=1,nufr(lm2l(lm),ias2is(ias))
                wnklm_exact(lm,ir,ias,ispn)=wnklm_exact(lm,ir,ias,ispn)+&
                  m_wann_unkmt(lm,io,ias,ispn,ik)*ufr(ir,lm2l(lm),io,ias2ic(ias))
              enddo !io
            enddo !lm
          enddo !ir
! convert to spherical coordinates
          call zgemm('T','N',mt_ntp,nrmt(ias2is(ias)),lmmaxvr,zone,mt_ylmf,&
            lmmaxvr,wnklm_exact(1,1,ias,ispn),lmmaxvr,zzero,&
            wnktp_exact(1,1,ias,ispn),mt_ntp)
        enddo !ias
      enddo !ispn
      wnkir_exact=zzero
      do ispn=1,nspinor
        do ig=1,m_ngknr(ik)
          wnkir_exact(igfft(m_igkignr(ig,ik)),ispn)=m_wann_unkit(ig,ispn,ik)/sqrt(omega)
        enddo
        call zfftifc(3,ngrid,1,wnkir_exact(1,ispn))
      enddo
      
      zt1=zzero
      zt2=zzero
      do ispn=1,nspinor
        zt1=zt1+s_zfinp(.true.,.false.,lmmaxvr,ngrtot,wnklm_exact(1,1,1,ispn),&
          wnklm_exact(1,1,1,ispn),wnkir_exact(1,ispn),wnkir_exact(1,ispn),zt2)
      enddo
      zprodlm_exact(1,j,ik)=zt1
      zprodlm_exact(2:3,j,ik)=zt2(1:2)
      
      zt1=zzero
      zt2=zzero
      do ispn=1,nspinor
        zt1=zt1+s_zfinp(.false.,.false.,mt_ntp,ngrtot,wnktp_exact(1,1,1,ispn),&
          wnktp_exact(1,1,1,ispn),wnkir_exact(1,ispn),wnkir_exact(1,ispn),zt2)
      enddo
      zprodtp_exact(1,j,ik)=zt1
      zprodtp_exact(2:3,j,ik)=zt2(1:2)
    
    enddo !ik
  enddo !j

  if (mpi_grid_root()) then
    open(210,file=trim(adjustl(fname)),form="formatted",status="replace")
    do ik=1,nkpt
      write(210,'(" ik : ",I4)')ik
      do j=1,sic_wantran%nwan
        n=sic_wantran%iwan(j)
        write(210,'("  n : ",I4,6X," exact lm <W_nk|W_nk> : ",3G18.10)')&
          n,dreal(zprodlm_exact(1,j,ik)),dreal(zprodlm_exact(2,j,ik)),dreal(zprodlm_exact(3,j,ik))
        write(210,'("  n : ",I4,6X," exact tp <W_nk|W_nk> : ",3G18.10)')&
          n,dreal(zprodtp_exact(1,j,ik)),dreal(zprodtp_exact(2,j,ik)),dreal(zprodtp_exact(3,j,ik))
      enddo
      write(210,*)
    enddo !ik
    close(210)
  endif
endif


if (itest.eq.2) then
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    call elk_load_wann_unk(n)
    do ik=1,nkpt
      wnktp_bt=zzero
      wnkir_bt=zzero
      do ias=1,natmtot
        is=ias2is(ias)
        ia=ias2ia(ias)
        do itloc=1,ntrloc
          it=mpi_grid_map(sic_blochsum%ntr,dim2,loc=itloc)
          expikt=exp(-zi*dot_product(vkc(:,ik),sic_blochsum%vtc(:,it)))
          do ir=1,nrmt(is)
            do itp=1,mt_ntp
              x(:)=mt_spx(:,itp)*spr(ir,is)+atposc(:,ia,is)+&
                   sic_blochsum%vtc(:,it)
              call ws_reduce(x)
              call s_get_wanval(x,wanval,rmax=sic_wan_rwsmax)
              wnktp_bt(itp,ir,ias,:)=wnktp_bt(itp,ir,ias,:)+expikt*wanval(:)
            enddo !itp
          enddo !ir
        enddo !itloc
      enddo !ias
      do it=1,sic_blochsum%ntr
        expikt=exp(-zi*dot_product(vkc(:,ik),sic_blochsum%vtc(:,it))) 
        do ir=1,ngrtot
          x(:)=vgrc(:,ir)+sic_blochsum%vtc(:,it)
          call ws_reduce(x)
          call s_get_wanval(x,wanval,rmax=sic_wan_rwsmax)
          wnkir_bt(ir,:)=wnkir_bt(ir,:)+expikt*wanval(:)
        enddo
      enddo
      zt1=zzero
      zt2=zzero
      do ispn=1,nspinor
        zt1=zt1+s_zfinp(.false.,.false.,mt_ntp,ngrtot,wnktp_bt(1,1,1,ispn),&
          wnktp_bt(1,1,1,ispn),wnkir_bt(1,ispn),wnkir_bt(1,ispn),zt2)
      enddo
      zprodtp_bt(1,j,ik)=zt1
      zprodtp_bt(2:3,j,ik)=zt2(1:2)

    enddo !ik
  enddo !j
  if (mpi_grid_root()) then
    open(210,file=trim(adjustl(fname)),form="formatted",status="replace")
    do ik=1,nkpt
      write(210,'(" ik : ",I4)')ik
      do j=1,sic_wantran%nwan
        n=sic_wantran%iwan(j)
        write(210,'("  n : ",I4,6X," bt tp    <W_nk|W_nk> : ",3G18.10)')&
          n,dreal(zprodtp_bt(1,j,ik)),dreal(zprodtp_bt(2,j,ik)),dreal(zprodtp_bt(3,j,ik))
        enddo
      write(210,*)
    enddo !ik
    close(210)
  endif
endif

if (itest.eq.3) then
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    do ik=1,nkpt
      wnktp_bt=zzero
      wnkir_bt=zzero
      do ias=1,natmtot
        is=ias2is(ias)
        ia=ias2ia(ias)
        do itloc=1,ntrloc
          it=mpi_grid_map(sic_blochsum%ntr,dim2,loc=itloc)
          expikt=exp(-zi*dot_product(vkc(:,ik),sic_blochsum%vtc(:,it)))
          do ir=1,nrmt(is)
            do itp=1,mt_ntp
              x(:)=mt_spx(:,itp)*spr(ir,is)+atposc(:,ia,is)+&
                   sic_blochsum%vtc(:,it)
              call ws_reduce(x)
              call s_spinor_func_val(x,s_wlm(1,1,1,j),wanval)
              wnktp_bt(itp,ir,ias,:)=wnktp_bt(itp,ir,ias,:)+expikt*wanval(:)
            enddo !itp
          enddo !ir
        enddo !itloc
      enddo !ias
      do it=1,sic_blochsum%ntr
        expikt=exp(-zi*dot_product(vkc(:,ik),sic_blochsum%vtc(:,it))) 
        do ir=1,ngrtot
          x(:)=vgrc(:,ir)+sic_blochsum%vtc(:,it)
          call ws_reduce(x)
          call s_spinor_func_val(x,s_wlm(1,1,1,j),wanval)
          wnkir_bt(ir,:)=wnkir_bt(ir,:)+expikt*wanval(:)
        enddo
      enddo
      zt1=zzero
      zt2=zzero
      do ispn=1,nspinor
        zt1=zt1+s_zfinp(.false.,.false.,mt_ntp,ngrtot,wnktp_bt(1,1,1,ispn),&
          wnktp_bt(1,1,1,ispn),wnkir_bt(1,ispn),wnkir_bt(1,ispn),zt2)
      enddo
      zprodtp_bt(1,j,ik)=zt1
      zprodtp_bt(2:3,j,ik)=zt2(1:2)

    enddo !ik
  enddo !j

  if (mpi_grid_root()) then
    open(210,file=trim(adjustl(fname)),form="formatted",status="replace")
    do ik=1,nkpt
      write(210,'(" ik : ",I4)')ik
      do j=1,sic_wantran%nwan
        n=sic_wantran%iwan(j)
        write(210,'("  n : ",I4,6X," bt tp    <W_nk|W_nk> : ",3G18.10)')&
          n,dreal(zprodtp_bt(1,j,ik)),dreal(zprodtp_bt(2,j,ik)),dreal(zprodtp_bt(3,j,ik))
        enddo
      write(210,*)
    enddo !ik
    close(210)
  endif
endif

deallocate(wnklm_exact)
deallocate(wnkir_exact)
deallocate(wnktp_exact)
deallocate(wnktp_bt)
deallocate(wnkir_bt)

deallocate(zprodlm_exact)
deallocate(zprodtp_exact)
deallocate(zprodtp_bt)


return
end
