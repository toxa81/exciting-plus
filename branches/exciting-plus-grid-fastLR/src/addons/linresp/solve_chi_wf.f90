subroutine solve_chi_wf(igq0,vcgq,w,nnzme,inzme,vcwan,f_response_)
use modmain
implicit none
integer, intent(in) :: igq0
real(8), intent(in) :: vcgq(ngvecchi)
complex(8), intent(in) :: w
integer, intent(in) :: nnzme
integer, intent(in) :: inzme(2,nnzme)
complex(8), intent(in) :: vcwan(nnzme,nnzme)
complex(8), intent(out) :: f_response_(nf_response)

complex(8), allocatable :: mtrx1(:,:)
complex(8), allocatable :: mtrx2(:,:)
complex(8), allocatable :: mtrx3(:,:)
complex zt1,zt2

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
    i3=itridxwan(i1,i2)
    if (i3.ne.-1) then
      mtrx1(i,j)=chi0wan(n1,n2,i3)
    endif
  enddo
enddo

mtrx3=zzero
do i1=1,nnzme
  mtrx3(i1,i1)=zone
enddo
call zgemm('N','N',nnzme,nnzme,nnzme,dcmplx(-1.d0,0.d0), &
  vcwan,nnzme,mtrx1,nnzme,dcmplx(1.d0,0.d0),mtrx3,nnzme)
call invzge(mtrx3,nnzme)
call zgemm('N','N',nnzme,nnzme,nnzme,dcmplx(1.d0,0.d0), &
  mtrx1,nnzme,mtrx3,nnzme,dcmplx(0.d0,0.d0),mtrx2,nnzme)

! zt1 is chi0wf
zt1=zzero
! zt2 is chiwf
zt2=zzero
do i=1,nnzme
  do j=1,nnzme
    i1=inzme(1,i)
    n1=inzme(2,i)
    i2=inzme(1,j)
    n2=inzme(2,j)
    i3=itridxwan(i1,i2)
    if (i3.ne.-1) then
      zt1=zt1+megqwan(n1,i1,igq0)*chi0wan(n1,n2,i3)*dconjg(megqwan(n2,i2,igq0))
    endif
    zt2=zt2+megqwan(n1,i1,igq0)*mtrx2(i,j)*dconjg(megqwan(n2,i2,igq0))
  enddo
enddo

deallocate(mtrx1)
deallocate(mtrx2)
deallocate(mtrx3)

f_response_(f_chi0_wann)=zt1
f_response_(f_chi_wann)=zt2

f_response_(f_epsilon_eff_wann)=1.d0/(1.d0+(vcgq(igq0)**2)*f_response_(f_chi_wann))
f_response_(f_sigma_wann)=zi*dreal(w)*(zone-f_response_(f_epsilon_eff_wann))/fourpi
f_response_(f_loss_wann)=1.d0/f_response_(f_epsilon_eff_wann)



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

