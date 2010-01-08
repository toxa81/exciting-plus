subroutine genwu(ng,chi0,vcgq,qnm)
implicit none
integer, intent(in) :: ng
complex(8), intent(in) :: chi0(ng,ng)
real(8), intent(in) :: vcgq(ng)
character(*), intent(in) :: qnm

complex(8), allocatable :: epsilon(:,:)
complex(8), allocatable :: mtrx1(:,:)
integer ig1,ig2,n1
character*100 fw

allocate(epsilon(ng,ng))
allocate(mtrx1(ng,ng))

!
!! rpa kernel
!mtrx1=dcmplx(0.d0,0.d0)
!do i=1,ng
!  mtrx1(i,i)=vcgq(i)**2
!enddo
!! compute matrix epsilon=1-chi0*v
!epsilon=dcmplx(0.d0,0.d0)
!do i=1,ng
!  epsilon(i,i)=dcmplx(1.d0,0.d0)
!enddo
!call zgemm('N','N',ng,ng,ng,dcmplx(-1.d0,0.d0),chi0,ng,mtrx1,ng,&
!  dcmplx(1.d0,0.d0),epsilon,ng)
!! invert epsilon matrix
!call invzge(epsilon,ng)
!! compute chi=epsilon^-1 * chi0
!call zgemm('N','N',ng,ng,ng,dcmplx(1.d0,0.d0),epsilon,ng,chi0,ng,&
!  dcmplx(0.d0,0.d0),mtrx1,ng)

! compute screened Coulomb potential using "symmetrized" dielectric function
do ig1=1,ng
  do ig2=1,ng
    epsilon(ig1,ig2)=-vcgq(ig1)*chi0(ig1,ig2)*vcgq(ig2)
  enddo
  epsilon(ig1,ig1)=dcmplx(1.d0,0.d0)+epsilon(ig1,ig1)
enddo
call invzge(epsilon,ng)
do ig1=1,ng
  do ig2=1,ng
    mtrx1(ig1,ig2)=vcgq(ig1)*epsilon(ig1,ig2)*vcgq(ig2)
  enddo
enddo


if (ng.gt.10) then
  n1=10
else
  n1=ng
endif
fw=trim(qnm)//"_w__.txt"
open(170,file=trim(fw),status='replace',form='formatted')
write(170,'("Screened W matrix")')
write(170,'("real part")')
do ig1=1,n1
  write(170,'(100F12.6)')(dreal(mtrx1(ig1,ig2)),ig2=1,n1)
enddo
write(170,'("imag part")')
do ig1=1,n1
  write(170,'(100F12.6)')(dimag(mtrx1(ig1,ig2)),ig2=1,n1)
enddo
close(170)

deallocate(epsilon,mtrx1)

return 
end