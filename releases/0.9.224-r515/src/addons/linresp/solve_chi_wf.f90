subroutine solve_chi_wf(ntr1,ntr2,itridx,nwfme,nnzme,inzme,mewf2, &
  !mewf4,mtrx_v,fourpiq0,epsilonwf,losswf,chi0wf,chiwf,igq0)
  mewf4,mtrx_v,fourpiq0,epsilonwf,chi0wf,chiwf,igq0)
use modmain
implicit none
integer, intent(in) :: ntr1
integer, intent(in) :: ntr2
integer, intent(in) :: itridx(ntr1,ntr1)
integer, intent(in) :: nwfme
integer, intent(in) :: nnzme
integer, intent(in) :: inzme(2,nnzme)
complex(8), intent(in) :: mewf2(nwfme,ntr1,ngvecme)
complex(8), intent(in) :: mewf4(nwfme,nwfme,ntr2)
complex(8), intent(in) :: mtrx_v(nnzme,nnzme)
real(8), intent(in) :: fourpiq0
complex(8), intent(out) :: epsilonwf
!complex(8), intent(out) :: losswf
complex(8), intent(out) :: chi0wf
complex(8), intent(out) :: chiwf
integer, intent(in) :: igq0

complex(8), allocatable :: mtrx1(:,:)
complex(8), allocatable :: mtrx2(:,:)
complex(8), allocatable :: mtrx3(:,:)
complex zt1

integer i1,i2,n1,n2,i,j,it2,i3,ig

complex(8), external :: zdotu

allocate(mtrx1(nnzme,nnzme))
allocate(mtrx2(nnzme,nnzme))
allocate(mtrx3(nnzme,nnzme))

mtrx1=zzero
do i=1,nnzme
  do j=1,nnzme
    i1=inzme(1,i)
    n1=inzme(2,i)
    i2=inzme(1,j)
    n2=inzme(2,j)
    i3=itridx(i1,i2)
    if (i3.ne.-1) then
      mtrx1(i,j)=mewf4(n1,n2,i3)
    endif
  enddo
enddo

mtrx3=zzero
do i1=1,nnzme
  mtrx3(i1,i1)=zone
enddo
call zgemm('N','N',nnzme,nnzme,nnzme,dcmplx(-1.d0,0.d0), &
  mtrx_v,nnzme,mtrx1,nnzme,dcmplx(1.d0,0.d0),mtrx3,nnzme)
call invzge(mtrx3,nnzme)
call zgemm('N','N',nnzme,nnzme,nnzme,dcmplx(1.d0,0.d0), &
  mtrx1,nnzme,mtrx3,nnzme,dcmplx(0.d0,0.d0),mtrx2,nnzme)

chi0wf=zzero
do i=1,nnzme
  do j=1,nnzme
    i1=inzme(1,i)
    n1=inzme(2,i)
    i2=inzme(1,j)
    n2=inzme(2,j)
    i3=itridx(i1,i2)
    if (i3.ne.-1) then
      chi0wf=chi0wf+mewf2(n1,i1,igq0)*mewf4(n1,n2,i3)*dconjg(mewf2(n2,i2,igq0))
    endif
    chiwf=chiwf+mewf2(n1,i1,igq0)*mtrx2(i,j)*dconjg(mewf2(n2,i2,igq0))
  enddo
enddo
! epsilon_eff
epsilonwf=1.d0/(1.d0+fourpiq0*chiwf)
! loss function
!losswf=1.d0/epsilonwf

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

