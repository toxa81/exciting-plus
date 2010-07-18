subroutine genuscrn(iq)
use modmain
use mod_addons_q
implicit none
integer, intent(in) :: iq
integer iwloc,nwloc,iw,n,n1,ig1,ig2,i
real(8) v2(3),vtc(3)
complex(8) zt1
complex(8), allocatable :: vscrn(:,:)
complex(8), allocatable :: megqwan1(:,:)
complex(8), allocatable :: expiqt(:)
complex(8), allocatable :: zm1(:,:)
complex(8), allocatable :: zm2(:,:)
allocate(vscrn(ngvecme,ngvecme))
allocate(megqwan1(ngvecme,nwann))                                                                                                     
allocate(expiqt(nmegqwan))
allocate(zm1(nwann,nwann))
allocate(zm2(nwann,ngvecme))

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
do iwloc=1,nwloc
  iw=mpi_grid_map(lr_nw,dim_k,loc=iwloc)
  call genvscrn(iq,chi0loc(1,1,iwloc),vscrn)
  call zgemm('C','N',nwann,ngvecme,ngvecme,zone,megqwan1,ngvecme,&
    vscrn,ngvecme,zzero,zm2,nwann)
  call zgemm('N','N',nwann,nwann,ngvecme,zone,zm2,nwann,megqwan1,ngvecme,&
    zzero,zm1,nwann)
  do i=1,nmegqwan
    n=imegqwan(1,i)
    n1=imegqwan(2,i)
!    v2=dble(imegqwan(3:5,i))
!    call r3mv(avec,v2,vtc)
!    zt1=zzero
!    do ig1=1,ngvecme
!      do ig2=1,ngvecme
!        zt1=zt1+dconjg(megqwan1(ig1,n))*vscrn(ig1,ig2)*&
!          megqwan1(ig2,n1)
!      enddo
!    enddo
!    uscrnwan(i,iwloc)=uscrnwan(i,iwloc)+zt1*exp(-zi*dot_product(vqc(:,iq),vtc))
    uscrnwan(i,iwloc)=uscrnwan(i,iwloc)+zm1(n,n1)*expiqt(i)
  enddo
enddo
deallocate(megqwan1)
deallocate(chi0loc)
deallocate(expiqt)
deallocate(zm1)
deallocate(zm2)
return
end