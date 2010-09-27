subroutine genuscrn(iq)
use modmain
use mod_addons_q
implicit none
integer, intent(in) :: iq
integer iwloc,nwloc,iw,n,n1,i,ig,nmegqwanloc,iloc
real(8) v2(3),vtc(3)
complex(8), allocatable :: vscrn(:,:)
complex(8), allocatable :: megqwan1(:,:)
complex(8), allocatable :: expiqt(:)
complex(8), allocatable :: zm1(:,:)
complex(8), allocatable :: zm2(:,:)
complex(8), allocatable :: krnl(:,:)
complex(8), allocatable :: epsilon(:,:)
complex(8), allocatable :: chi(:,:)

call papi_timer_start(pt_uscrn)

allocate(vscrn(ngvecme,ngvecme))
allocate(megqwan1(ngvecme,nwann))                                                                                                     
allocate(expiqt(nmegqwan))
allocate(zm1(nwann,nwann))
allocate(zm2(nwann,ngvecme))
allocate(krnl(ngvecme,ngvecme))
allocate(epsilon(ngvecme,ngvecme))
allocate(chi(ngvecme,ngvecme))
krnl=zzero
do ig=1,ngvecme
  krnl(ig,ig)=vhgq(ig,iq)
enddo
do i=1,nmegqwan
  v2=dble(imegqwan(3:5,i))
  call r3mv(avec,v2,vtc)
  expiqt(i)=exp(-zi*dot_product(vqc(:,iq),vtc))
enddo
do n=1,nwann
  megqwan1(:,n)=megqwan(idxmegqwan(n,n,0,0,0),:)
enddo
! distribute frequency points over 1-st dimension
nwloc=mpi_grid_map(lr_nw,dim_k)
! distribute Wannier transitions over 3-rd dimension
nmegqwanloc=mpi_grid_map(nmegqwan,dim_b)
do iwloc=1,nwloc
  iw=mpi_grid_map(lr_nw,dim_k,loc=iwloc)
! broadcast chi0
  call mpi_grid_bcast(chi0loc(1,1,iwloc),ngvecme*ngvecme,dims=(/dim_b/))
  call genvscrn(iq,iw,chi0loc(1,1,iwloc),krnl,vscrn,epsilon,chi)
  call zgemm('C','N',nwann,ngvecme,ngvecme,zone,megqwan1,ngvecme,&
    vscrn,ngvecme,zzero,zm2,nwann)
  call zgemm('N','N',nwann,nwann,ngvecme,zone,zm2,nwann,megqwan1,ngvecme,&
    zzero,zm1,nwann)
  do iloc=1,nmegqwanloc
    i=mpi_grid_map(nmegqwan,dim_b,loc=iloc)
    n=imegqwan(1,i)
    n1=imegqwan(2,i)
    uscrnwan(i,iwloc)=uscrnwan(i,iwloc)+zm1(n,n1)*expiqt(i)
  enddo
enddo
deallocate(megqwan1)
deallocate(chi0loc)
deallocate(expiqt)
deallocate(zm1)
deallocate(zm2)
deallocate(krnl,epsilon,chi)

call papi_timer_stop(pt_uscrn)

return
end