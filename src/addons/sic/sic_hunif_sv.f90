subroutine sic_hunif_sv(ikloc,evecfd)
use modmain
use mod_sic
use mod_seceqn
implicit none
integer, intent(in) :: ikloc
complex(8), intent(inout) :: evecfd(nspinor*nmatmax,nstsv)
!
integer ik,ig,ist,ispn,lm,ir,is,ic,ias,j,io,l,j1,j2
complex(8) zt1
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfsvmt(:,:,:,:,:)
complex(8), allocatable :: wfsvit(:,:,:)
!
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)   
allocate(wfsvmt(lmmaxapw,nufrmax,natmtot,nspinor,nstsv))
allocate(wfsvit(ngkmax,nspinor,nstsv))
allocate(apwalm(ngkmax,lmmaxapw,apwordmax,natmtot))
wfsvmt=zzero
wfsvit=zzero
call genapwalm(ngk(1,ik),gkc(:,1,ikloc),tpgkc(:,:,1,ikloc),&
  sfacgk(:,:,1,ikloc),apwalm)
call genwfsvc(lmaxapw,lmmaxapw,ngk(1,ik),nstsv,apwalm,evecfd,wfsvmt,wfsvit) 
deallocate(apwalm)
sic_wb(:,:,:,ikloc)=zzero
do j=1,sic_wantran%nwan
  do ispn=1,nspinor
    do ias=1,natmtot
      is=ias2is(ias)
      ic=ias2ic(ias)
      do lm=1,lmmaxapw
        l=lm2l(lm)
        do io=1,nufr(l,is)
          zt1=zzero
          do ir=1,nrmt(is)
            zt1=zt1+dconjg(s_wkmtlm(ir,lm,ias,ispn,j,ikloc))*ufr(ir,l,io,ic)*mt_rw(ir,is)
          enddo
          do ist=1,nstsv
            sic_wb(j,ist,1,ikloc)=sic_wb(j,ist,1,ikloc)+zt1*wfsvmt(lm,io,ias,ispn,ist)
          enddo
        enddo
      enddo
    enddo
    do ist=1,nstsv
      do ig=1,ngk(1,ik)
        sic_wb(j,ist,1,ikloc)=sic_wb(j,ist,1,ikloc)+&
          dconjg(s_wkit(ig,ispn,j,ikloc))*wfsvit(ig,ispn,ist)
      enddo
    enddo
  enddo
enddo

write(*,*)"in sic_hunif_sv"
do j1=1,sic_wantran%nwan
  do j2=1,sic_wantran%nwan
    zt1=zzero
    do j=1,nstsv
      zt1=zt1+dconjg(sic_wb(j1,j,1,ikloc))*sic_wb(j2,j,1,ikloc)
    enddo
    write(*,*)j1,j2,zt1
  enddo
enddo
deallocate(wfsvmt,wfsvit)
call pstop

return
end subroutine
