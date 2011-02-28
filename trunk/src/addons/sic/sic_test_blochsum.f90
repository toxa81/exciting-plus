subroutine sic_test_blochsum(itest,twan,fname)
use modmain
use mod_sic
implicit none
integer, intent(in) :: itest
logical, intent(in) :: twan
character*(*), intent(in) :: fname
integer ik,ikloc,j,n,jas,it,ias,is,ia,ir,itp,ispn,ntrloc,itloc,lm
integer ig
real(8) x(3),t1
real(8), allocatable :: tp(:,:)
complex(8), allocatable :: zprod(:,:,:)
complex(8) zt1,zt2(2),expikt,wanval(nspinor)
complex(8), allocatable :: wgk(:,:)
real(8), allocatable :: jl(:,:)
complex(8),allocatable :: ylmgk(:)
complex(8), allocatable :: wankmt(:,:,:,:)
complex(8), allocatable :: wankir(:,:)
complex(8), allocatable :: wantp(:,:,:,:)
complex(8), allocatable :: wanir(:,:)
!
allocate(zprod(3,sic_wantran%nwan,nkpt))
allocate(wgk(ngkmax,nspinor))
allocate(jl(s_nr,0:lmaxwan))
allocate(ylmgk(lmmaxwan))
allocate(wankmt(mt_ntp,nrmtmax,natmtot,nspinor))
allocate(wankir(ngrtot,nspinor))
allocate(wantp(mt_ntp,nrmtmax,natmtot,nspinor))
allocate(wanir(ngrtot,nspinor))


zprod=zzero
ntrloc=mpi_grid_map(sic_orbitals%ntr,dim2)
! test Bloch sum using exact values of WFs
if (itest.eq.1) then
  do ik=1,1 !nkpt
    do j=1,sic_wantran%nwan
      n=sic_wantran%iwan(j)
      wankmt=zzero
      wankir=zzero
      do ias=1,natmtot
        is=ias2is(ias)
        ia=ias2ia(ias)
        do itloc=1,ntrloc
          it=mpi_grid_map(sic_orbitals%ntr,dim2,loc=itloc)
          expikt=exp(-zi*dot_product(vkc(:,ik),sic_orbitals%vtc(:,it)))
          wantp=zzero
          do ir=1,nrmt(is)
            do itp=1,mt_ntp
              x(:)=mt_spx(:,itp)*spr(ir,is)+atposc(:,ia,is)+&
                   sic_orbitals%vtc(:,it)
              call s_get_wanval(.true.,n,x,wanval)
              wantp(itp,ir,ias,:)=wanval(:)
            enddo !itp
          enddo !ir
          call mpi_grid_reduce(wantp(1,1,1,1),mt_ntp*nrmtmax*natmtot*nspinor,&
            dims=(/dim_k/))
          wankmt(:,:,:,:)=wankmt(:,:,:,:)+expikt*wantp(:,:,:,:)
        enddo !itloc
      enddo !ias
      call mpi_grid_reduce(wankmt(1,1,1,1),mt_ntp*nrmtmax*natmtot*nspinor,&
        dims=(/dim2/))
      do itloc=1,ntrloc
        it=mpi_grid_map(sic_orbitals%ntr,dim2,loc=itloc)
        expikt=exp(-zi*dot_product(vkc(:,ik),sic_orbitals%vtc(:,it)))
        wanir=zzero
        do ir=1,ngrtot
          x(:)=vgrc(:,ir)+sic_orbitals%vtc(:,it)
          call s_get_wanval(twan,n,x,wanval)
          wanir(ir,:)=wanval(:)
        enddo
        call mpi_grid_reduce(wanir(1,1),ngrtot*nspinor,dims=(/dim_k/))
        wankir(:,:)=wankir(:,:)+expikt*wanir(:,:)
      enddo !itloc
      call mpi_grid_reduce(wankir(1,1),ngrtot*nspinor,dims=(/dim2/),&
        all=.true.)
      zt1=zzero
      zt2=zzero
      do ispn=1,nspinor
        zt1=zt1+s_zfinp(.false.,mt_ntp,wankmt(1,1,1,ispn),wankmt(1,1,1,ispn),&
          wankir(1,ispn),wankir(1,ispn),zt2)
      enddo
      zprod(1,j,ik)=zt1
      zprod(2:3,j,ik)=zt2(1:2)
    enddo !j
  enddo !ik
endif
if (itest.eq.2) then
  do ikloc=1,1 !nkptloc
    ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
    do j=1,sic_wantran%nwan
      n=sic_wantran%iwan(j)
      wgk=zzero
      do ig=1,ngk(1,ik)
! generate Bessel functions j_l(|G+k|x)
        do ir=1,s_nr
          call sbessel(lmaxwan,gkc(ig,1,ikloc)*s_r(ir),jl(ir,0))
        enddo
        call genylm(lmaxwan,tpgkc(1,ig,1,ikloc),ylmgk)
        do ispn=1,nspinor
          do ir=1,s_nr
            do lm=1,lmmaxwan
              wgk(ig,ispn)=wgk(ig,ispn)+fourpi*(zi**lm2l(lm))*jl(ir,lm2l(lm))*&
                dconjg(s_pwanlm(lm,ir,ispn,j))*dconjg(ylmgk(lm))*s_rw(ir)/sqrt(omega)
            enddo
          enddo
        enddo
      enddo !ig
      wankmt=zzero
      wankir=zzero
      do ias=1,natmtot
        is=ias2is(ias)
        ia=ias2ia(ias)
        do itloc=1,ntrloc
          it=mpi_grid_map(sic_orbitals%ntr,dim2,loc=itloc)
          expikt=exp(-zi*dot_product(vkc(:,ik),sic_orbitals%vtc(:,it)))
          do ir=1,nrmt(is)
            do itp=1,mt_ntp
              x(:)=mt_spx(:,itp)*spr(ir,is)+atposc(:,ia,is)+&
                   sic_orbitals%vtc(:,it)-wanpos(:,n)
              do ispn=1,nspinor
                wankmt(itp,ir,ias,ispn)=wankmt(itp,ir,ias,ispn)+&
                    expikt*s_func_val(x,s_wanlm(1,1,ispn,j))
              enddo
            enddo !itp
          enddo !ir
        enddo !itloc
      enddo !ias
      call mpi_grid_reduce(wankmt(1,1,1,1),mt_ntp*nrmtmax*natmtot*nspinor,&
        dims=(/dim2/),all=.true.)
      do ispn=1,nspinor
        do ig=1,ngk(1,ik)
          wankir(igfft(igkig(ig,1,ikloc)),ispn)=wgk(ig,ispn)/sqrt(omega)
        enddo
        call zfftifc(3,ngrid,1,wankir(1,ispn))
      enddo
      zt1=zzero
      zt2=zzero
      do ispn=1,nspinor
        zt1=zt1+s_zfinp(.false.,mt_ntp,wankmt(1,1,1,ispn),wankmt(1,1,1,ispn),&
          wankir(1,ispn),wankir(1,ispn),zt2)
      enddo
      zprod(1,j,ik)=zt1
      zprod(2:3,j,ik)=zt2(1:2)
    enddo !j
  enddo !ikloc 
endif
call mpi_grid_reduce(zprod(1,1,1),3*sic_wantran%nwan*nkpt,dims=(/dim_k/))
if (mpi_grid_root()) then
  open(210,file=trim(adjustl(fname)),form="formatted",status="replace")
  do ik=1,nkpt
    write(210,'(" ik : ",I4)')ik
    do j=1,sic_wantran%nwan
      n=sic_wantran%iwan(j)
      write(210,'("  n : ",I4,6X," <W_nk|W_nk> : ",3G18.10)')&
        n,dreal(zprod(1,j,ik)),dreal(zprod(2,j,ik)),dreal(zprod(3,j,ik))
    enddo
    write(210,*)
  enddo !ik
  close(210)
endif
deallocate(zprod)
deallocate(wgk,jl,ylmgk)
deallocate(wankmt,wankir)
deallocate(wantp,wanir)
return
end
