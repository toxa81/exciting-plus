subroutine sic_blochsum_mt
use modmain
use mod_sic
implicit none
integer ik,ikloc,j,n,it,ias,is,ia,ir,itp,ispn,ntrloc,itloc
real(8) x(3)
complex(8) zt1(nspinor),zt2(nspinor)
complex(8), allocatable :: expikt(:,:)
!
s_wankmt=zzero
s_wvkmt=zzero
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
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  do ias=1,natmtot
    is=ias2is(ias)
    ia=ias2ia(ias)
    do itloc=1,ntrloc
      it=mpi_grid_map(sic_blochsum%ntr,dim2,loc=itloc)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(itp,ispn,ikloc,x,zt1,zt2)
      do ir=1,nrmt(is)
        do itp=1,mt_ntp
          x(:)=mt_spx(:,itp)*spr(ir,is)+atposc(:,ia,is)+&
               sic_blochsum%vtc(:,it)-wanpos(:,n)
          call s_func_val2(x,s_wanlm(1,1,1,j),s_wvlm(1,1,1,j),&
            zt1,zt2)
          do ispn=1,nspinor
            do ikloc=1,nkptloc
              s_wankmt(itp,ir,ias,ispn,j,ikloc)=&
                s_wankmt(itp,ir,ias,ispn,j,ikloc)+&
                expikt(ikloc,itloc)*zt1(ispn)
              s_wvkmt(itp,ir,ias,ispn,j,ikloc)=&
                s_wvkmt(itp,ir,ias,ispn,j,ikloc)+&
                expikt(ikloc,itloc)*zt2(ispn)
            enddo !ikloc
          enddo !ispn
        enddo !itp
      enddo !ir
!$OMP END PARALLEL DO
    enddo !itloc
  enddo !ias
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
call timer_stop(90)
if (mpi_grid_root()) then
  write(*,'("[sic_blochsum_mt] total time : ",F12.4," sec.")')timer_get_value(90)
endif
return
end
