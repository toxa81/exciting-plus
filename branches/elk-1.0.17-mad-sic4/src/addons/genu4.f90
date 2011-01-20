subroutine genu4(iq,nwloc,ntrloc)
use modmain
use mod_addons_q
use mod_wannier
use mod_expigqr
use mod_linresp
implicit none
integer, intent(in) :: iq
integer, intent(in) :: nwloc
integer, intent(in) :: ntrloc
integer iwloc,iw,n,n1,i,ig,itloc,vtl(3),j,it
real(8) v2(3),vtc(3)
complex(8), allocatable :: vscrn(:,:)
complex(8), allocatable :: megqwan2(:,:)
complex(8), allocatable :: megqwan3(:,:)
complex(8), allocatable :: zm1(:,:)
complex(8), allocatable :: zm2(:,:)
complex(8), allocatable :: krnl(:,:)
complex(8), allocatable :: epsilon(:,:)
complex(8), allocatable :: chi(:,:)
complex(8) zt1

if (screenu4) call genchi0(iq)

call papi_timer_start(pt_uscrn)
allocate(vscrn(ngvecme,ngvecme))
allocate(krnl(ngvecme,ngvecme))
allocate(epsilon(ngvecme,ngvecme))
allocate(chi(ngvecme,ngvecme))
allocate(zm1(megqwantran%nwt,ngvecme))
allocate(zm2(megqwantran%nwt,megqwantran%nwt))
krnl=zzero
do ig=1,ngvecme
  krnl(ig,ig)=vhgq(ig,iq)
enddo
allocate(megqwan2(ngvecme,megqwantran%nwt))   
allocate(megqwan3(ngvecme,megqwantran%nwt))   
! compute megqwan2=<W_n|e^{+i(G+q)x}|W_n'T'> and also rearrange megqwan
do i=1,megqwantran%nwt
  n=megqwantran%iwt(1,i)
  n1=megqwantran%iwt(2,i)
  vtl(:)=megqwantran%iwt(3:5,i)
  v2=dble(vtl)
  call r3mv(avec,v2,vtc)
  zt1=exp(-zi*dot_product(vqc(:,iq),vtc))
  j=megqwantran%iwtidx(n1,n,-vtl(1),-vtl(2),-vtl(3))
  if (j.le.0) then
    write(*,'("Error(genu4) wrong index of matrix element")')
    write(*,'(" n,n1,vtl : ",5I4)')n,n1,vtl
    call pstop
  endif
  megqwan2(:,i)=dconjg(megqwan(j,:))*zt1
  megqwan3(:,i)=megqwan(i,:)
enddo
! compute 4-index U
! TODO: comments with formulas
do iwloc=1,nwloc
  iw=mpi_grid_map(lr_nw,dim_k,loc=iwloc)
! broadcast chi0
  if (screenu4) then
    call mpi_grid_bcast(chi0loc(1,1,iwloc),ngvecme*ngvecme,dims=(/dim_b/))
    call genvscrn(iq,chi0loc(1,1,iwloc),krnl,vscrn,epsilon,chi)
  else
    vscrn=krnl
  endif
  call zgemm('T','N',megqwantran%nwt,ngvecme,ngvecme,zone,megqwan2,ngvecme,&
    vscrn,ngvecme,zzero,zm1,megqwantran%nwt)
  call zgemm('N','N',megqwantran%nwt,megqwantran%nwt,ngvecme,zone,zm1,&
    megqwantran%nwt,megqwan3,ngvecme,zzero,zm2,megqwantran%nwt)
  do itloc=1,ntrloc
    it=mpi_grid_map(megqwantran%ntr,dim_b,loc=itloc)
    v2=dble(megqwantran%vtr(:,it))
    call r3mv(avec,v2,vtc)
    zt1=exp(-zi*dot_product(vqc(:,iq),vtc))/omega/nkptnr
    call zaxpy((megqwantran%nwt)**2,zt1,zm2(1,1),1,u4(1,1,itloc,iwloc),1)
  enddo
enddo
deallocate(megqwan2)
deallocate(megqwan3)
deallocate(zm1)
deallocate(zm2)
deallocate(krnl,epsilon,chi)
call papi_timer_stop(pt_uscrn)
return
end