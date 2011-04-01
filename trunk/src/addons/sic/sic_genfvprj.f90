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
integer ias,ist,j,j1
complex(8) zt2(2)
complex(8), allocatable :: wfmt_(:,:)
complex(8), allocatable :: wfmt(:,:,:)
complex(8), allocatable :: wb(:,:,:,:)
complex(8), allocatable :: om(:,:)
!
sic_wb(:,:,:,ikloc)=zzero
sic_wvb(:,:,:,ikloc)=zzero
if (.not.tsic_wv) return 
if (debug_level.ge.4) then
  allocate(wb(3,sic_wantran%nwan,nstfv,nspinor))
  wb=zzero
  allocate(om(sic_wantran%nwan,sic_wantran%nwan))
endif
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
      zt2=zzero
      sic_wb(j,ist,ispn,ikloc)=s_zfinp(.false.,.true.,mt_ntp,ngk(1,ik),&
        s_wankmt(1,1,1,ispn,j,ikloc),wfmt,s_wankir(1,ispn,j,ikloc),&
        evecfv(1,ist,1),zt2)
      if (debug_level.ge.4) then  
        wb(1,j,ist,ispn)=sic_wb(j,ist,ispn,ikloc)
        wb(2:3,j,ist,ispn)=zt2(:)
      endif
      sic_wvb(j,ist,ispn,ikloc)=s_zfinp(.false.,.true.,mt_ntp,ngk(1,ik),&
        s_wvkmt(1,1,1,ispn,j,ikloc),wfmt,s_wvkir(1,ispn,j,ikloc),&
        evecfv(1,ist,1))
    enddo
  enddo
enddo !ist
!$OMP END DO
deallocate(wfmt_,wfmt)
!$OMP END PARALLEL
call timer_stop(t_sic_genfvprj)
if (debug_level.ge.4) then
  call dbg_open_file
  write(fdbgout,*)
  write(fdbgout,'("[sic_genfvprj] ik : ",I4,"    vkl : ",3F10.6)')ik,vkl(:,ik)
  do ist=1,nstfv
    write(fdbgout,'("  ist : ",I4)')ist
    do j=1,sic_wantran%nwan 
      write(fdbgout,'("    j  : ",I4)')j
      do ispn=1,nspinor
        write(fdbgout,'("     total : ",F12.8," (",2F12.8,")")')abs(wb(1,j,ist,ispn)),&
          dreal(wb(1,j,ist,ispn)),dimag(wb(1,j,ist,ispn))
        write(fdbgout,'("        mt : ",F12.8," (",2F12.8,")")')abs(wb(2,j,ist,ispn)),&
          dreal(wb(2,j,ist,ispn)),dimag(wb(2,j,ist,ispn))
        write(fdbgout,'("        it : ",F12.8," (",2F12.8,")")')abs(wb(3,j,ist,ispn)),&
          dreal(wb(3,j,ist,ispn)),dimag(wb(3,j,ist,ispn))
      enddo
    enddo
  enddo
  om=zzero
  do j=1,sic_wantran%nwan
    do j1=1,sic_wantran%nwan
      do ispn=1,nspinor
        do ist=1,nstfv
          om(j,j1)=om(j,j1)+wb(1,j,ist,ispn)*dconjg(wb(1,j1,ist,ispn))
        enddo
      enddo
    enddo
  enddo
  write(fdbgout,'("  overlap matrix")')
  do j=1,sic_wantran%nwan
    write(fdbgout,'(255F12.6)')(abs(om(j,j1)),j1=1,sic_wantran%nwan)
  enddo
  call dbg_close_file
  deallocate(wb,om)
endif
return
end
