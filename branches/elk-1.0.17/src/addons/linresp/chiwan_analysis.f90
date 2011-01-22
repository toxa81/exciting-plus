subroutine chiwan_analysis(fout,chi0wan,vcwan,iq,iw,w,thresh)
use modmain
use mod_addons_q

implicit none
integer, intent(in) :: fout
complex(8), intent(in) :: chi0wan(nmegqwan,nmegqwan)
complex(8), intent(in) :: vcwan(nmegqwan,nmegqwan)
integer, intent(in) :: iq
integer, intent(in) :: iw
real(8), intent(in) :: w
real(8), intent(in) :: thresh

complex(8), allocatable :: mtrx1(:,:)
complex(8), allocatable :: mtrx2(:,:)
complex(8), allocatable ::  zt2(:)
real(8), allocatable :: rt2(:)
integer, allocatable :: lindx_i(:)
integer, allocatable :: lindx_j(:)
character*100 :: str,c1
integer nmegqwan2
integer i,j,li,lj
real(8) t1
integer it1, it2

allocate(mtrx1(nmegqwan,nmegqwan))
allocate(mtrx2(nmegqwan,nmegqwan))

! square
nmegqwan2=nmegqwan*nmegqwan

allocate(lindx_i(nmegqwan2))
allocate(lindx_j(nmegqwan2))
allocate(zt2(nmegqwan2))
allocate(rt2(nmegqwan2))
lindx_i=0
lindx_j=0
zt2=zzero
rt2=0.0

str="omega index: "
write(c1,'(i6)')iw
str=trim(str)//adjustl(c1)
write(fout,'(A)')str
str="omega value (eV): "
write(c1,'(f7.3)')w
str=trim(str)//adjustl(c1)
write(fout,'(A)')str
write(fout,'(A)')""

! commemt:
! compute chi0_GqGq using the Wannier functions expansion
!  chi0=\sum_{T,T'}\sum_{n,m,n',m'} A^{*}_{nmT}(q,Gq)*chi0wan(n,m,n',m',T-T')*A_{n'm'T'}(q,Gq)
!  where A_{nmT}(q,G)=<n,0|e^{-i(G+q)x}|m,T>

mtrx1=zzero
do i=1,nmegqwan
  mtrx1(i,i)=zone
enddo
! mtrx1 = 1-v*chi0
call zgemm('N','N',nmegqwan,nmegqwan,nmegqwan,-zone,vcwan,nmegqwan,chi0wan,&
  nmegqwan,zone,mtrx1,nmegqwan)
call invzge(mtrx1,nmegqwan)
! mtrx2 = chi0*(1-v*chi0)^-1
call zgemm('N','N',nmegqwan,nmegqwan,nmegqwan,zone,chi0wan,nmegqwan,mtrx1,&
  nmegqwan,zzero,mtrx2,nmegqwan)
! zt1 is chi0_GqGq_wan
!zt1=zzero
! zt2 is chi_GqGq_wan
zt2=zzero
li=0
do i=1,nmegqwan
  do j=1,nmegqwan
!    zt1=zt1+megqwan(i,iig0q)*chi0wan(i,j)*dconjg(megqwan(j,iig0q))
    li=li+1
    lindx_i(li)=i   
    lindx_j(li)=j   
    zt2(li)=megqwan(i,iig0q)*mtrx2(i,j)*dconjg(megqwan(j,iig0q))
    rt2(li)=abs(zt2(li))
!    if (abs(zt2).gt.thresh) then
!      write(fout,'(i6i6f10.4e20.8)')i,j,w,abs(zt2)
!    endif
  enddo
enddo

! sort by abs(zt1)
do li=1,nmegqwan2
  do lj=1,li
    if (rt2(lj).lt.rt2(li)) then
      t1=rt2(li)
      rt2(li)=rt2(lj)
      rt2(lj)=t1
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
    if (rt2(li).gt.0.000001) then
      i=lindx_i(li)
      j=lindx_j(li)
      write(fout,'(i6i6d15.4d15.4d15.4d20.8)')i,j,abs(megqwan(i,iig0q)),abs(dconjg(megqwan(j,iig0q))),abs(chi0wan(i,j)),rt2(li)
    endif
enddo



deallocate(mtrx1)
deallocate(mtrx2)
return
end

