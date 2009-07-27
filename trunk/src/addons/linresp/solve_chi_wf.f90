subroutine solve_chi_wf(ntr1,ntr2,itridx,nwfme,ndim,mewf2,mewf4,mtrx_v,chi0wf,chiwf)
use modmain
implicit none
integer, intent(in) :: ntr1
integer, intent(in) :: ntr2
integer, intent(in) :: itridx(ntr1,ntr1)
integer, intent(in) :: nwfme
integer, intent(in) :: ndim
complex(8), intent(in) :: mewf2(nwfme,ntr1,ngvecme)
complex(8), intent(in) :: mewf4(nwfme,nwfme,ntr2)
complex(8), intent(in) :: mtrx_v(ndim,ndim)
complex(8), intent(out) :: chi0wf
complex(8), intent(out) :: chiwf


complex(8), allocatable :: mtrx1(:,:)
complex(8), allocatable :: mtrx2(:,:)
complex(8), allocatable :: mtrx3(:,:)
complex zt1

integer i1,i2,n1,n2,i,j,it2,i3,ig

complex(8), external :: zdotu

allocate(mtrx1(ndim,ndim))
allocate(mtrx2(ndim,ndim))
allocate(mtrx3(ndim,ndim))

mtrx1=zzero
do i1=1,ntr1
  do i2=1,ntr1
    i3=itridx(i1,i2)
    if (i3.ne.-1) then
      do n1=1,nwfme
        do n2=1,nwfme
          mtrx1((i1-1)*nwfme+n1,(i2-1)*nwfme+n2)=mewf4(n1,n2,i3)
        enddo
      enddo
    endif
  enddo !i2
enddo !i1

mtrx3=zzero
do i1=1,ndim
  mtrx3(i1,i1)=zone
enddo
call zgemm('N','N',ndim,ndim,ndim,dcmplx(-1.d0,0.d0), &
  mtrx_v,ndim,mtrx1,ndim,dcmplx(1.d0,0.d0),mtrx3,ndim)
call invzge(mtrx3,ndim)
call zgemm('N','N',ndim,ndim,ndim,dcmplx(1.d0,0.d0), &
  mtrx1,ndim,mtrx3,ndim,dcmplx(0.d0,0.d0),mtrx2,ndim)

chi0wf=zzero
do it2=1,ntr2
  do i=1,ntr1
    do j=1,ntr1
      if (itridx(i,j).eq.it2) then
        do n2=1,nwfme
          zt1=zdotu(nwfme,mewf4(1,n2,it2),1,mewf2(1,i,1),1)
          chi0wf=chi0wf+zt1*dconjg(mewf2(n2,j,1))
        enddo
      endif
    enddo
  enddo
enddo !it2

chiwf=zzero
do i=1,ntr1
  do j=1,ntr1
    do n2=1,nwfme
      zt1=zdotu(nwfme,mtrx2((i-1)*nwfme+1,(j-1)*nwfme+n2),1,mewf2(1,i,1),1)
      chiwf=chiwf+zt1*dconjg(mewf2(n2,j,1))
    enddo
  enddo
enddo

deallocate(mtrx1)
deallocate(mtrx2)
deallocate(mtrx3)

!  if (lwannopt) then
!    do it1=1,ntr1
!      do it2=1,ntr1
!        do n1=1,nwann*nwann
!          do n2=1,nwann*nwann
!            epswf(ie)=epswf(ie)+mewf4(n1,n2,itridx(it1,it2))*mewfx(1,n1,it1)*dconjg(mewfx(1,n2,it2))
!          enddo
!        enddo
!      enddo
!    enddo 
!  endif


return
end

