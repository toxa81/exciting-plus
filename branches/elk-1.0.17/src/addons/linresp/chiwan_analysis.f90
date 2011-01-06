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
complex zt1,zt2
character*100 :: str,c1

integer i,j

allocate(mtrx1(nmegqwan,nmegqwan))
allocate(mtrx2(nmegqwan,nmegqwan))

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
zt1=zzero
! zt2 is chi_GqGq_wan
zt2=zzero
do i=1,nmegqwan
  do j=1,nmegqwan
!    zt1=zt1+megqwan(i,iig0q)*chi0wan(i,j)*dconjg(megqwan(j,iig0q))
    zt2=zt2+megqwan(i,iig0q)*mtrx2(i,j)*dconjg(megqwan(j,iig0q))
    if (abs(zt2).gt.thresh) then
      write(fout,'(i6i6f10.4e20.8)')i,j,w,abs(zt2)
    endif
  enddo
enddo
deallocate(mtrx1)
deallocate(mtrx2)
return
end

