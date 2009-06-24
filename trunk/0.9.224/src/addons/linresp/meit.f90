subroutine zrhoftit(nme0,ime0,ngknr1,ngknr2,igkignr1,igkignr2, &
  igkq,wfsvit1,wfsvit2,zrhofc0,gvit)
use modmain
implicit none
! arguments
integer, intent(in) :: nme0
integer, intent(in) :: ime0(3,nmemax)
integer, intent(in) :: ngknr1
integer, intent(in) :: ngknr2
integer, intent(in) :: igkq
integer, intent(in) :: igkignr1(ngkmax)
integer, intent(in) :: igkignr2(ngkmax)
complex(8), intent(in) :: wfsvit1(ngkmax,nspinor,nstsv)
complex(8), intent(in) :: wfsvit2(ngkmax,nspinor,nstsv)
complex(8), intent(inout) :: zrhofc0(ngvecme,nmemax)
complex(8), intent(in) :: gvit(intgv(1,1):intgv(1,2),intgv(2,1):intgv(2,2),intgv(3,1):intgv(3,2))

complex(8), allocatable :: mit(:,:)
complex(8), allocatable :: a(:,:,:) 
complex(8), allocatable :: zrhofc_tmp(:,:)

integer ig,ist1,ist2,i,ispn,ispn2,j,ist
integer idx0,bs,idx_g1,idx_g2,igp,ifg,ir,i1,i2
complex(8) zt1
complex(8), external :: zdotu,zdotc

allocate(zrhofc_tmp(ngvecme,nmemax))
zrhofc_tmp=zzero

allocate(a(ngknr1,nstsv,nspinor))
allocate(mit(ngknr1,ngknr2))

call idxbos(ngvecme,mpi_dims(2),mpi_x(2)+1,idx0,bs)
idx_g1=idx0+1
idx_g2=idx0+bs

do ig=idx_g1,idx_g2
  call genpwit2(ngknr1,ngknr2,igkignr1,igkignr2,-(ivg(:,ig+gvecme1-1)+ivg(:,igkq)),gvit,mit)
  a=zzero
  do ispn=1,nspinor
    do i=1,nstsv
      call zgemv('N',ngknr1,ngknr2,zone,mit,ngknr1,wfsvit2(1,ispn,i),1,zzero,a(1,i,ispn),1)
    enddo
  enddo
  do ispn=1,nspinor
    if (lrtype.eq.0) then
      ispn2=ispn
    endif
    if (lrtype.eq.1) then
      ispn2=3-ispn
    endif
    do i=1,nme0
      ist1=ime0(1,i)
      ist2=ime0(2,i)
      zrhofc_tmp(ig,i)=zrhofc_tmp(ig,i)+&
        zdotc(ngknr1,wfsvit1(1,ispn,ist1),1,a(1,ist2,ispn2),1)
    enddo
  enddo
enddo !ig  
deallocate(mit,a)

if (mpi_dims(2).gt.1) then
  call d_reduce_cart(comm_cart_010,.false.,zrhofc_tmp,2*ngvecme*nmemax)
endif
zrhofc0=zrhofc0+zrhofc_tmp
deallocate(zrhofc_tmp)

return
end