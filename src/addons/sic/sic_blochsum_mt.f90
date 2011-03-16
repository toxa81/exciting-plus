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
ntrloc=mpi_grid_map(sic_orbitals%ntr,dim2)
allocate(expikt(nkptloc,ntrloc))
do itloc=1,ntrloc
  it=mpi_grid_map(sic_orbitals%ntr,dim2,loc=itloc)
  do ikloc=1,nkptloc
    ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
    expikt(ikloc,itloc)=&
      exp(-zi*dot_product(vkc(:,ik),sic_orbitals%vtc(:,it)))
  enddo
enddo
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  do ias=1,natmtot
    is=ias2is(ias)
    ia=ias2ia(ias)
    do itloc=1,ntrloc
      it=mpi_grid_map(sic_orbitals%ntr,dim2,loc=itloc)
      do ir=1,nrmt(is)
        do itp=1,mt_ntp
          x(:)=mt_spx(:,itp)*spr(ir,is)+atposc(:,ia,is)+&
               sic_orbitals%vtc(:,it)-wanpos(:,n)
          call s_func_val2(x,s_wanlm(1,1,1,j),s_wvlm(1,1,1,j),&
            zt1,zt2)
          do ispn=1,nspinor
            !zt1=s_func_val(x,s_wanlm(1,1,ispn,j))
            !zt2=s_func_val(x,s_wvlm(1,1,ispn,j))
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

!ntrloc=mpi_grid_map(sic_orbitals%ntr,dim2)
!do ikloc=1,nkptloc
!  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
!! make Bloch sums
!  do j=1,sic_wantran%nwan
!    n=sic_wantran%iwan(j)
!    do ias=1,natmtot
!      is=ias2is(ias)
!      ia=ias2ia(ias)
!      do itloc=1,ntrloc
!        it=mpi_grid_map(sic_orbitals%ntr,dim2,loc=itloc)
!        expikt=exp(-zi*dot_product(vkc(:,ik),sic_orbitals%vtc(:,it)))
!        do ir=1,nrmt(is)
!          do itp=1,mt_ntp
!            x(:)=mt_spx(:,itp)*spr(ir,is)+atposc(:,ia,is)+&
!                 sic_orbitals%vtc(:,it)-wanpos(:,n)
!            do ispn=1,nspinor
!              s_wankmt(itp,ir,ias,ispn,j,ikloc)=s_wankmt(itp,ir,ias,ispn,j,ikloc)+&
!                expikt*s_func_val(x,s_wanlm(1,1,ispn,j))
!            enddo
!          enddo
!        enddo
!      enddo !itloc
!    enddo !ias
!    call mpi_grid_reduce(s_wankmt(1,1,1,1,j,ikloc),mt_ntp*nrmtmax*natmtot*nspinor,&
!      dims=(/dim2/),all=.true.)
!  enddo !j
!enddo !ikloc
!if (mpi_grid_root()) then
!  open(210,file="blochsum_from_bt_wannier.dat",form="formatted",&
!    status="replace")
!  do ir=1,nrmt(1)
!    !t1=0.d0
!    !do itp=1,mt_ntp
!    !  t1=t1+abs(abs(s_wankmt(itp,ir,1,1,1,ikloc))-abs(wfmt(itp,ir,1)))
!    !enddo
!    itp=3
!    write(210,'(3G18.10)')spr(ir,1),dreal(s_wankmt(itp,ir,1,1,1,1)),&
!      dimag(s_wankmt(itp,ir,1,1,1,1))
!  enddo
!  close(210)
!endif
call timer_stop(90)
if (mpi_grid_root()) then
  write(*,'("[sic_blochsum_mt] total time : ",F12.4," sec.")')timer_get_value(90)
endif
return
end
