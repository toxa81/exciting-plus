subroutine genwu(iw,chi0,vcgq,qnm,vscr)
use modmain
implicit none
integer, intent(in) :: iw
complex(8), intent(in) :: chi0(ngvecme,ngvecme)
real(8), intent(in) :: vcgq(ngvecme)
character(*), intent(in) :: qnm
complex(8), intent(out) :: vscr(ngvecme,ngvecme)

complex(8), allocatable :: epsilon(:,:)
!complex(8), allocatable :: mtrx1(:,:)
integer ig1,ig2,n1,n2
integer i,j
character*100 fname
integer, external :: hash
character*12 c12

allocate(epsilon(ngvecme,ngvecme))
!allocate(mtrx1(ngvecme,ngvecme))

!! rpa kernel
!mtrx1=dcmplx(0.d0,0.d0)
!do i=1,ngvecme
!  mtrx1(i,i)=vcgq(i)**2
!enddo
!! compute matrix epsilon=1-chi0*v
!epsilon=dcmplx(0.d0,0.d0)
!do i=1,ngvecme
!  epsilon(i,i)=dcmplx(1.d0,0.d0)
!enddo
!call zgemm('N','N',ngvecme,ngvecme,ngvecme,dcmplx(-1.d0,0.d0),chi0,ngvecme,mtrx1,ngvecme,&
!  dcmplx(1.d0,0.d0),epsilon,ngvecme)
!! invert epsilon matrix
!call invzge(epsilon,ngvecme)
!! compute chi=epsilon^-1 * chi0
!call zgemm('N','N',ngvecme,ngvecme,ngvecme,dcmplx(1.d0,0.d0),epsilon,ngvecme,chi0,ngvecme,&
!  dcmplx(0.d0,0.d0),mtrx1,ngvecme)

! compute screened Coulomb potential using "symmetrized" dielectric function
do ig1=1,ngvecme
  do ig2=1,ngvecme
    epsilon(ig1,ig2)=-vcgq(ig1)*chi0(ig1,ig2)*vcgq(ig2)
  enddo
  epsilon(ig1,ig1)=zone+epsilon(ig1,ig1)
enddo
call invzge(epsilon,ngvecme)
do ig1=1,ngvecme
  do ig2=1,ngvecme
    vscr(ig1,ig2)=vcgq(ig1)*epsilon(ig1,ig2)*vcgq(ig2)
  enddo
enddo


! compute screened u
do ig1=1,ngvecme
  do ig2=1,ngvecme
    do n1=1,nwann
      do n2=1,nwann
        uscrnwan(n1,n2,iw)=uscrnwan(n1,n2,iw)+&
          dconjg(megqwan(idxmegqwan(n1,n1,0,0,0),ig1))*vscr(ig1,ig2)*&
          megqwan(idxmegqwan(n2,n2,0,0,0),ig2)
      enddo
    enddo
  enddo
enddo
! compute bare u
!ubare=zzero
if (iw.eq.1) then
  do ig1=1,ngvecme
    do n1=1,nwann
      do n2=1,nwann
        ubarewan(n1,n2)=ubarewan(n1,n2)+&
          dconjg(megqwan(idxmegqwan(n1,n1,0,0,0),ig1))*(vcgq(ig1)**2)*&
          megqwan(idxmegqwan(n2,n2,0,0,0),ig1)
      enddo
    enddo
  enddo
endif
deallocate(epsilon)
return 
end