subroutine sic_e0
use modmain
use modldapu
use mod_mpi_grid
use mod_sic
implicit none
complex(8), allocatable :: ene0(:,:,:,:,:)
integer n1,n2,ias,lm1,lm2,ispn1,ispn2,ikloc,ik,ispn,n,j1,j2

! note: we need H_{nn'}(k) to compute energies of WFs: E_n=<W_n|H|W_n>
!  which is computed as \sum_{k}H_{nn}(k)
if (.not.tsic_wv) then
  do ikloc=1,nkptloc
    ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
    call genwann_h(.false.,evalsv(1,ik),wann_c(1,1,ikloc),wann_h(1,1,ik),&
      wann_e(1,ik))
    do j1=1,sic_wantran%nwan
      n1=sic_wantran%iwan(j1)
      do j2=1,sic_wantran%nwan
        n2=sic_wantran%iwan(j1)
        sic_wann_h0k(j1,j2,ikloc)=wann_h(n1,n2,ik)
      enddo
    enddo
  enddo !ikloc
endif
allocate(ene0(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
ene0=zzero
! collinear (!!!) case only
do j1=1,sic_wantran%nwan
  n1=sic_wantran%iwan(j1)
   do j2=1,sic_wantran%nwan
    n2=sic_wantran%iwan(j2)
    if (wan_info(1,n1).eq.wan_info(1,n2)) then
      ias=wan_info(1,n1)
      lm1=wan_info(2,n1)
      lm2=wan_info(2,n2)
      ispn1=wan_info(3,n1)
      ispn2=wan_info(3,n2)
      do ikloc=1,nkptloc
        ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)        
        ene0(lm1,lm2,ispn1,ispn2,ias)=ene0(lm1,lm2,ispn1,ispn2,ias)+&
          sic_wann_h0k(j1,j2,ikloc)*wkpt(ik)
      enddo
    endif
  enddo
enddo 
call mpi_grid_reduce(ene0(1,1,1,1,1),&
  lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot,dims=(/dim_k/),all=.true.)
! convert from Rlm to Ylm basis
do ias=1,natmtot
  do ispn=1,nspinor
    call unimtrxt(lmmaxlu,yrlm_lps(1,1,ias),ene0(1,1,ispn,ispn,ias))
  enddo
enddo
! symmetrise matrix
call symdmat(lmaxlu,lmmaxlu,ene0)
! convert back to Rlm
do ias=1,natmtot
  do ispn=1,nspinor
    call unimtrxt(lmmaxlu,rylm_lps(1,1,ias),ene0(1,1,ispn,ispn,ias))
  enddo
enddo
do n=1,nwantot
  ias=wan_info(1,n)
  lm1=wan_info(2,n)
  ispn1=wan_info(3,n)
  sic_wann_e0(n)=ene0(lm1,lm1,ispn1,ispn1,ias)
enddo
if (mpi_grid_root()) then
  open(170,file="SIC_WANN_E0.OUT",form="FORMATTED",status="REPLACE")
  do n=1,nwantot
    write(170,'(G18.10)')sic_wann_e0(n)
  enddo
  close(170)
endif
deallocate(ene0)
return
end

