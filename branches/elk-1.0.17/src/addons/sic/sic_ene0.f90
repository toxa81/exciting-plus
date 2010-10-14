subroutine sic_e0
use modmain
use modldapu
use mod_mpi_grid
implicit none
complex(8), allocatable :: ene0(:,:,:,:,:)
integer n1,n2,ias,lm1,lm2,ispn1,ispn2,ikloc,ik,ispn,n

allocate(ene0(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))
ene0=zzero
! collinear (!!!) case only
do n1=1,nwann
  do n2=1,nwann
    if (iwann(1,n1).eq.iwann(1,n2)) then
      ias=iwann(1,n1)
      lm1=iwann(2,n1)
      lm2=iwann(2,n2)
      ispn1=iwann(3,n1)
      ispn2=iwann(3,n2)
      do ikloc=1,nkptloc
        ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)        
        ene0(lm1,lm2,ispn1,ispn2,ias)=ene0(lm1,lm2,ispn1,ispn2,ias)+&
          sic_wann_h0k(n1,n2,ikloc)*wkpt(ik)
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
do n=1,nwann
  ias=iwann(1,n)
  lm1=iwann(2,n)
  ispn1=iwann(3,n)
  sic_wann_e0(n)=ene0(lm1,lm1,ispn1,ispn1,ias)
enddo
if (mpi_grid_root()) then
  open(170,file="SIC_WANN_E0.OUT",form="FORMATTED",status="REPLACE")
  do n=1,nwann
    write(170,'(G18.10)')sic_wann_e0(n)
  enddo
  close(170)
endif
deallocate(ene0)
return
end

