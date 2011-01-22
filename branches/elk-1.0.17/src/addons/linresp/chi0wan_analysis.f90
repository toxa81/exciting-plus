subroutine chi0wan_analysis(fout,chi0wan,iq,iw,w,thresh)
use modmain
use mod_addons_q

! arguments
integer, intent(in) :: fout
complex(8), intent(in) :: chi0wan(nmegqwan,nmegqwan)
integer, intent(in) :: iq
integer, intent(in) :: iw
real(8), intent(in) :: w
real(8), intent(in) :: thresh
! local variables
integer i,j,li,nmegqwan2
integer, allocatable :: lindx_i(:)
integer, allocatable :: lindx_j(:)
complex(8), allocatable :: zt1(:)
real(8), allocatable :: rt1(:)
character*100 :: str,c1
real(8) t1
integer it1, it2

! square
nmegqwan2=nmegqwan*nmegqwan

allocate(lindx_i(nmegqwan2))
allocate(lindx_j(nmegqwan2))
allocate(zt1(nmegqwan2))
allocate(rt1(nmegqwan2))
lindx=0
lindx_i=0
lindx_j=0
zt1=zzero
rt1=0.0

str="omega index: "
write(c1,'(i6)')iw
str=trim(str)//adjustl(c1)
write(fout,'(A)')str
str="omega value (eV): "
write(c1,'(f7.3)')w
str=trim(str)//adjustl(c1)
write(fout,'(A)')""
write(fout,'(A)')str
str="iproc : "
write(c1,'(i5)')iproc
str=trim(str)//adjustl(c1)
write(fout,'(A)')str
write(fout,'(A)')""

! build array with a linear indx of A_{\lambda} 
! \chi_{\lambda \lambda'} A_{\lambda'}
li=0
do i=1,nmegqwan
  do j=1,nmegqwan
    ! linear index
    !li=N*(i-1)+j
    li=li+1
    lindx_i(li)=i   
    lindx_j(li)=j   
    zt1(li)=megqwan(i,iig0q)*chi0wan(i,j)*dconjg(megqwan(j,iig0q))
    rt1(li)=abs(zt1(li))
  enddo
enddo

!absolute value
!rt1=abs(zt1)

! sort by abs(zt1)
do li=1,nmegqwan2
  do lj=1,li
    if (rt1(lj).lt.rt1(li)) then
      t1=rt1(li)
      rt1(li)=rt1(lj)
      rt1(lj)=t1
      it1=lindx_i(li)
      it2=lindx_j(li)
      lindx_i(li)=lindx_i(lj)
      lindx_j(li)=lindx_j(lj)
      lindx_i(lj)=it1
      lindx_j(lj)=it2
    endif
  enddo
enddo

! output to screen/file
do li=1,nmegqwan2
    if (rt1(li).gt.0.000001) then
      i=lindx_i(li)
      j=lindx_j(li)
      write(fout,'(i6i6d15.4d15.4d15.4d20.8)')i,j,abs(megqwan(i,iig0q)),abs(dconjg(megqwan(j,iig0q))),abs(chi0wan(i,j)),rt1(li)
      !write(fout,'(i6i6e20.8e20.8e20.8e20.8)')i,j,abs(megqwan(i,iig0q)),abs(dconjg(megqwan(j,iig0q))),abs(chi0wan(i,j)),rt1(li)
    endif
enddo

deallocate(lindx_i)
deallocate(lindx_j)
deallocate(zt1)
deallocate(rt1)

end
