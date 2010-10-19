subroutine genu4(iq)
use modmain
use mod_addons_q
implicit none
integer, intent(in) :: iq
integer iwloc,nwloc,iw,n,n1,i,ig,ig1,itloc,iloc,vtl(3),j,it,ntmegqwanloc
real(8) v2(3),vtc(3)
complex(8), allocatable :: vscrn(:,:)
complex(8), allocatable :: megqwan2(:,:)
complex(8), allocatable :: megqwan3(:,:)
complex(8), allocatable :: expiqt(:)
complex(8), allocatable :: zm1(:,:)
complex(8), allocatable :: zm2(:,:)
complex(8), allocatable :: krnl(:,:)
complex(8), allocatable :: epsilon(:,:)
complex(8), allocatable :: chi(:,:)
complex(8) zt1

call papi_timer_start(pt_uscrn)

allocate(vscrn(ngvecme,ngvecme))
allocate(krnl(ngvecme,ngvecme))
allocate(epsilon(ngvecme,ngvecme))
allocate(chi(ngvecme,ngvecme))
allocate(zm1(nmegqwan,ngvecme))
allocate(zm2(nmegqwan,nmegqwan))
krnl=zzero
do ig=1,ngvecme
  krnl(ig,ig)=vhgq(ig,iq)
enddo
! distribute frequency points over 1-st dimension
nwloc=mpi_grid_map(lr_nw,dim_k)
! distribute translations along 3-nd dimention
ntmegqwanloc=mpi_grid_map(ntmegqwan,dim_b)
allocate(megqwan2(ngvecme,nmegqwan))   
allocate(megqwan3(ngvecme,nmegqwan))   
! compute megqwan2=<W_n|e^{+i(G+q)x}|W_n'T'> and also rearrange megqwan
do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  vtl(:)=imegqwan(3:5,i)
  v2=dble(vtl)
  call r3mv(avec,v2,vtc)
  zt1=exp(-zi*dot_product(vqc(:,iq),vtc))
  j=idxmegqwan(n1,n,-vtl(1),-vtl(2),-vtl(3))
  if (j.le.0) then
    write(*,'("Error(genuscrn) wrong index of matrix element")')
    write(*,'(" n,n1,vtl : ",5I4)')n,n1,vtl
    call pstop
  endif
  megqwan2(:,i)=dconjg(megqwan(j,:))*zt1
  megqwan3(:,i)=megqwan(i,:)
enddo
! compute 4-index U
do iwloc=1,nwloc
  iw=mpi_grid_map(lr_nw,dim_k,loc=iwloc)
! broadcast chi0
  call mpi_grid_bcast(chi0loc(1,1,iwloc),ngvecme*ngvecme,dims=(/dim_b/))
  call genvscrn(iq,iw,chi0loc(1,1,iwloc),krnl,vscrn,epsilon,chi)
  call zgemm('T','N',nmegqwan,ngvecme,ngvecme,zone,megqwan2,ngvecme,&
    vscrn,ngvecme,zzero,zm1,nmegqwan)
  call zgemm('N','N',nmegqwan,nmegqwan,ngvecme,zone,zm1,nmegqwan,megqwan3,&
    ngvecme,zzero,zm2,nmegqwan)
  do itloc=1,ntmegqwanloc
    it=mpi_grid_map(ntmegqwan,dim_b,loc=itloc)
    v2=dble(itmegqwan(:,it))
    call r3mv(avec,v2,vtc)
    zt1=exp(-zi*dot_product(vqc(:,iq),vtc))/omega/nkptnr
    call zaxpy(nmegqwan*nmegqwan,zt1,zm2(1,1),1,u4(1,1,itloc,iwloc),1)
  enddo
enddo

deallocate(megqwan2)
deallocate(megqwan3)
deallocate(chi0loc)
deallocate(zm1)
deallocate(zm2)
deallocate(krnl,epsilon,chi)

call papi_timer_stop(pt_uscrn)

return
end