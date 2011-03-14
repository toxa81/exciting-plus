subroutine sic_blochsum_mt
use modmain
use mod_sic
implicit none
integer ik,ikloc,j,n,it,ias,is,ia,ir,itp,ispn,ntrloc,itloc
integer itp1
real(8) x(3)
complex(8) expikt,zt1
!
s_wankmt=zzero
s_wvkmt=zzero
if (.not.tsic_wv) return
ntrloc=mpi_grid_map(sic_orbitals%ntr,dim2)
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
! make Bloch sums
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
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
              s_wankmt(itp,ir,ias,ispn,j,ikloc)=s_wankmt(itp,ir,ias,ispn,j,ikloc)+&
                expikt*s_func_val(x,s_wanlm(1,1,ispn,j))
            enddo
          enddo
        enddo
      enddo !itloc
    enddo !ias
    call mpi_grid_reduce(s_wankmt(1,1,1,1,j,ikloc),mt_ntp*nrmtmax*natmtot*nspinor,&
      dims=(/dim2/),all=.true.)
  enddo !j
enddo !ikloc
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
return
end
