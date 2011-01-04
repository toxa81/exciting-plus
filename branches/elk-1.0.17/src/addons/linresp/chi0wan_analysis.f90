subroutine chi0wan_analysis(fout,chi0wan,iq,w,thresh)
use modmain
use mod_addons_q

! arguments
integer, intent(in) :: fout
complex(8), intent(in) :: chi0wan(nmegqwan,nmegqwan)
integer, intent(in) :: iq
real(8), intent(in) :: w
real(8), intent(in) :: thresh
! local variables
integer i,j,li,nmegqwan2
integer, allocatable :: lindx(:)
integer, allocatable :: lindx_i(:)
integer, allocatable :: lindx_j(:)
complex(8), allocatable :: zt1(:)
real(8), allocatable :: rt1(:)

! square
nmegqwan2=nmegqwan*nmegqwan

allocate(lindx(nmegqwan2))
allocate(lindx_i(nmegqwan2))
allocate(lindx_j(nmegqwan2))
allocate(zt1(nmegqwan2))
allocate(rt1(nmegqwan2))
lindx=0
lindx_i=0
lindx_j=0
zt1=zzero
rt1=0.0

! build array with a linear indx of A_{\lambda} 
! \chi_{\lambda \lambda'} A_{\lambda'}
do i=1,nmegqwan
  do j=1,nmegqwan
    ! linear index
    li=N*(i-1)+j
    lindx_i(li)=i   
    lindx_j(li)=j   
    zt1(li)=megqwan(i,iig0q)*chi0wan(i,j)*dconjg(megqwan(j,iig0q))
    if (abs(zt1(li)).gt.thresh) then
      write(fout,'(i6i6f10.4e20.8)')i,j,w,abs(zt1(li))
!      write(*,'(i6i6f10.4e20.8)')i,j,w,abs(zt1(li))
    endif
  enddo
enddo

!absolute value
rt1=abs(zt1)

! sort by abs(zt1)
!call sortidx(nmegqwan2,rt1,lindx)

! output to screen/file
!do i=1,nmegqwan2
!  j=lindx(i)
!  write(*,'(i6i6e20.8)')lindx_i(j),lindx_j(j),rt1(j)
!enddo

deallocate(lindx)
deallocate(lindx_i)
deallocate(lindx_j)
deallocate(zt1)
deallocate(rt1)

end
