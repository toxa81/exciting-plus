subroutine sic_blochsum_mt
use modmain
use mod_sic
use mod_ws
implicit none
integer ik,ikloc,j,n,it,ias,is,ia,ir,itp,ispn,ntrloc,itloc,ig
real(8) x(3),t1(3)
complex(8) zt1(nspinor),zt2(nspinor),z1
complex(8), allocatable :: expikt(:,:)
complex(8), allocatable :: zfir(:,:,:)
integer, external :: hash
!
write(*,*)"hash(wanlm)=",hash(s_wanlm,lmmaxwan*s_nr*nspinor*sic_wantran%nwan*16)
s_wankmt=zzero
s_wvkmt=zzero
s_wankir=zzero
if (.not.tsic_wv) return
call timer_start(90,reset=.true.)
ntrloc=mpi_grid_map(sic_blochsum%ntr,dim2)
allocate(expikt(nkptloc,ntrloc))
do itloc=1,ntrloc
  it=mpi_grid_map(sic_blochsum%ntr,dim2,loc=itloc)
  do ikloc=1,nkptloc
    ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
    expikt(ikloc,itloc)=&
      exp(-zi*dot_product(vkc(:,ik),sic_blochsum%vtc(:,it)))
  enddo
enddo
allocate(zfir(ngrtot,nspinor,nkptloc))
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  do ias=1,natmtot
    is=ias2is(ias)
    ia=ias2ia(ias)
    do itloc=1,ntrloc
      it=mpi_grid_map(sic_blochsum%ntr,dim2,loc=itloc)
!!!!!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(itp,ispn,ikloc,x,zt1,zt2)
      do ir=1,nrmt(is)
        do itp=1,mt_ntp
          x(:)=mt_spx(:,itp)*spr(ir,is)+atposc(:,ia,is)+&
               sic_blochsum%vtc(:,it)-wanpos(:,n)
          call ws_reduce(x,sic_wan_rmax,t1)
          !t1=t1(:)-sic_blochsum%vtc(:,it)
          call s_func_val2(x,s_wanlm(1,1,1,j),s_wvlm(1,1,1,j),&
            zt1,zt2,sic_wan_rmax)
          do ispn=1,nspinor
            do ikloc=1,nkptloc
              ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
              s_wankmt(itp,ir,ias,ispn,j,ikloc)=&
                s_wankmt(itp,ir,ias,ispn,j,ikloc)+&
                exp(-zi*dot_product(vkc(:,ik),t1(:)))*zt1(ispn)
              s_wvkmt(itp,ir,ias,ispn,j,ikloc)=&
                s_wvkmt(itp,ir,ias,ispn,j,ikloc)+&
                expikt(ikloc,itloc)*zt2(ispn)
            enddo !ikloc
          enddo !ispn
        enddo !itp
      enddo !ir
!!!!!$OMP END PARALLEL DO
    enddo !itloc
  enddo !ias
  zfir=zzero
  do itloc=1,ntrloc
    it=mpi_grid_map(sic_blochsum%ntr,dim2,loc=itloc)
    do ir=1,ngrtot
      x(:)=vgrc(:,ir)+sic_blochsum%vtc(:,it)-wanpos(:,n)
      call ws_reduce(x,sic_wan_rmax,t1)
      !t1=t1(:)-sic_blochsum%vtc(:,it)
      call s_func_val2(x,s_wanlm(1,1,1,j),s_wvlm(1,1,1,j),&
          zt1,zt2,sic_wan_rmax)
      do ispn=1,nspinor
        do ikloc=1,nkptloc
          ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
          zfir(ir,ispn,ikloc)=zfir(ir,ispn,ikloc)+&
              exp(-zi*dot_product(vkc(:,ik),t1(:)))*zt1(ispn)
        enddo !ikloc
      enddo !ispn
    enddo !ir
  enddo
  do ikloc=1,nkptloc
    z1=zzero
    do ispn=1,nspinor
      do ir=1,ngrtot
        z1=z1+zfir(ir,ispn,ikloc)*dconjg(zfir(ir,ispn,ikloc))*cfunir(ir)*omega/dble(ngrtot)  
      enddo
    enddo
    write(*,*)"ikloc=",ikloc,"zt1=",z1
  enddo
  do ispn=1,nspinor
    do ikloc=1,nkptloc
      ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
      do ig=1,ngk(1,ik)
        z1=zzero
        do ir=1,ngrtot
          z1=z1+exp(zi*dot_product(vgkc(:,ig,1,ikloc),vgrc(:,ir)))*zfir(ir,ispn,ikloc)*&
            cfunir(ir)*omega/dble(ngrtot)/sqrt(omega)
        enddo
        s_wankir(ig,ispn,j,ikloc)=z1
      enddo
      !call zfftifc(3,ngrid,-1,zfir(1,ispn,ikloc))
      !zt1=zzero
      !do ig=1,ngk(1,ik)
      !  s_wankir(ig,ispn,j,ikloc)=zfir(igfft(igkig(ig,1,ikloc)),ispn,ikloc)
      !  zt1=zt1+s_wankir(ig,ispn,j,ikloc)*dconjg(s_wankir(ig,ispn,j,ikloc))
      !enddo
      !write(*,*)"ikloc=",ikloc,"zt1=",zt1
    enddo
  enddo
enddo !j
do ikloc=1,nkptloc
  do j=1,sic_wantran%nwan
    call mpi_grid_reduce(s_wankmt(1,1,1,1,j,ikloc),&
      mt_ntp*nrmtmax*natmtot*nspinor,dims=(/dim2/),all=.true.)
    call mpi_grid_reduce(s_wvkmt(1,1,1,1,j,ikloc),&
      mt_ntp*nrmtmax*natmtot*nspinor,dims=(/dim2/),all=.true.)
  enddo !j
enddo !ikloc
deallocate(expikt)
deallocate(zfir)
call timer_stop(90)
if (mpi_grid_root()) then
  write(*,'("[sic_blochsum_mt] total time : ",F12.4," sec.")')timer_get_value(90)
endif
return
end
