subroutine solve_chi_wan(igq0,vcgq,w,nmegqwan2,imegqwan2,megqwan2,vcwan,chi0wan,f_response_)
use modmain
implicit none
integer, intent(in) :: igq0
real(8), intent(in) :: vcgq(ngvecme)
complex(8), intent(in) :: w
integer, intent(in) :: nmegqwan2
integer, intent(in) :: imegqwan2(2,nmegqwan2)
complex(8), intent(in) :: megqwan2(nmegqwan,ntrmegqwan,ngvecme)
complex(8), intent(in) :: vcwan(nmegqwan2,nmegqwan2)
complex(8), intent(in) :: chi0wan(nmegqwan,nmegqwan,ntrchi0wan)
complex(8), intent(out) :: f_response_(nf_response)

complex(8), allocatable :: mtrx1(:,:)
complex(8), allocatable :: mtrx2(:,:)
complex(8), allocatable :: mtrx3(:,:)
complex zt1,zt2

integer i1,i2,n1,n2,i,j,it2,i3,ig

complex(8), external :: zdotu

allocate(mtrx1(nmegqwan2,nmegqwan2))
allocate(mtrx2(nmegqwan2,nmegqwan2))
allocate(mtrx3(nmegqwan2,nmegqwan2))

mtrx1=zzero
do i=1,nmegqwan2
  do j=1,nmegqwan2
    i1=imegqwan2(1,i)
    n1=imegqwan2(2,i)
    i2=imegqwan2(1,j)
    n2=imegqwan2(2,j)
    i3=itridxwan(i1,i2)
    mtrx1(i,j)=chi0wan(n1,n2,i3)
  enddo
enddo

mtrx3=zzero
do i1=1,nmegqwan2
  mtrx3(i1,i1)=zone
enddo
call zgemm('N','N',nmegqwan2,nmegqwan2,nmegqwan2,dcmplx(-1.d0,0.d0), &
  vcwan,nmegqwan2,mtrx1,nmegqwan2,dcmplx(1.d0,0.d0),mtrx3,nmegqwan2)
call invzge(mtrx3,nmegqwan2)
call zgemm('N','N',nmegqwan2,nmegqwan2,nmegqwan2,dcmplx(1.d0,0.d0), &
  mtrx1,nmegqwan2,mtrx3,nmegqwan2,dcmplx(0.d0,0.d0),mtrx2,nmegqwan2)

! zt1 is chi0wf
zt1=zzero
! zt2 is chiwf
zt2=zzero
do i=1,nmegqwan2
  do j=1,nmegqwan2
    i1=imegqwan2(1,i)
    n1=imegqwan2(2,i)
    i2=imegqwan2(1,j)
    n2=imegqwan2(2,j)
    i3=itridxwan(i1,i2)
    zt1=zt1+megqwan2(n1,i1,igq0)*chi0wan(n1,n2,i3)*dconjg(megqwan2(n2,i2,igq0))
    zt2=zt2+megqwan2(n1,i1,igq0)*mtrx2(i,j)*dconjg(megqwan2(n2,i2,igq0))
  enddo
enddo
f_response_(f_chi0_wann)=zt1
f_response_(f_chi_wann)=zt2
f_response_(f_epsilon_eff_wann)=1.d0/(1.d0+(vcgq(igq0)**2)*f_response_(f_chi_wann))
f_response_(f_sigma_wann)=zi*dreal(w)*(zone-f_response_(f_epsilon_eff_wann))/fourpi
f_response_(f_loss_wann)=1.d0/f_response_(f_epsilon_eff_wann)
deallocate(mtrx1)
deallocate(mtrx2)
deallocate(mtrx3)
return
end

