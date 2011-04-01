subroutine sic_genfvprj(ikloc,evecfv,apwalm)
use modmain
use mod_sic
implicit none
! arguments
integer, intent(in) :: ikloc
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
! local variables
integer ik,ispn
integer ias,ist,j
complex(8) zt2(2)
complex(8), allocatable :: wfmt_(:,:)
complex(8), allocatable :: wfmt(:,:,:)
!
sic_wb(:,:,:,ikloc)=zzero
sic_wvb(:,:,:,ikloc)=zzero
if (.not.tsic_wv) return 
call timer_start(t_sic_genfvprj)
!do ikloc=1,nkptloc
!  do j=1,sic_wantran%nwan
!    n=sic_wantran%iwan(j)
!    do ispn=1,nspinor
!      do ist=1,nstfv
!        istsv=ist+(ispn-1)*nstfv
!        do i=1,nstsv
!          sic_wb_tmp(j,ist,ispn,ikloc)=sic_wb_tmp(j,ist,ispn,ikloc)+&
!            dconjg(wann_c(n,i,ikloc)*evecsvloc(istsv,i,ikloc))
!        enddo !j
!      enddo !ispn
!    enddo !i
!  enddo !n
!enddo !ikloc  
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(wfmt_,wfmt,ias,j,ispn)
allocate(wfmt_(lmmaxvr,nrmtmax))
allocate(wfmt(mt_ntp,nrmtmax,natmtot))
!$OMP DO
do ist=1,nstfv
  wfmt=zzero
! generate first-variational wave function
  do ias=1,natmtot
    call wavefmt(1,lmaxvr,ias2is(ias),ias2ia(ias),ngk(1,ik),apwalm,&
      evecfv(1,ist,1),lmmaxvr,wfmt_)
! convert to spherical coordinates
    call zgemm('T','N',mt_ntp,nrmt(ias2is(ias)),lmmaxvr,zone,mt_ylmf,&
      lmmaxvr,wfmt_,lmmaxvr,zzero,wfmt(1,1,ias),mt_ntp)
  enddo
  do j=1,sic_wantran%nwan
    do ispn=1,nspinor
      sic_wb(j,ist,ispn,ikloc)=s_zfinp(.false.,.true.,mt_ntp,ngk(1,ik),&
        s_wankmt(1,1,1,ispn,j,ikloc),wfmt,s_wankir(1,ispn,j,ikloc),&
        evecfv(1,ist,1))
      sic_wvb(j,ist,ispn,ikloc)=s_zfinp(.false.,.true.,mt_ntp,ngk(1,ik),&
        s_wvkmt(1,1,1,ispn,j,ikloc),wfmt,s_wvkir(1,ispn,j,ikloc),&
        evecfv(1,ist,1))
    enddo
  enddo
  !write(*,'("    i : ",I4,"  re,im,abs : ",3F10.6," ! mt,it : ",&
  !        2F10.6,",",2F10.6)')ist,dreal(sic_wb(1,ist,1,ik)),dimag(sic_wb(1,ist,1,ik)),&
  !        abs(sic_wb(1,ist,1,ik)),dreal(zt2(1)),dimag(zt2(1)),&
  !        dreal(zt2(2)),dimag(zt2(2))
enddo !ist
!$OMP END DO
!write(*,*)
deallocate(wfmt_,wfmt)
!$OMP END PARALLEL
call timer_stop(t_sic_genfvprj)
return
end