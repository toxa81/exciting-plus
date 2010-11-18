subroutine solve_chi_wan(iq,w,vcwan,chi0wan,f_response_)
use modmain
use mod_addons_q
implicit none
integer, intent(in) :: iq
complex(8), intent(in) :: w
complex(8), intent(in) :: vcwan(nmegqwan,nmegqwan)
complex(8), intent(in) :: chi0wan(nmegqwan,nmegqwan)
complex(8), intent(out) :: f_response_(nf_response)

complex(8), allocatable :: mtrx1(:,:)
complex(8), allocatable :: mtrx2(:,:)
complex zt1,zt2

integer i,j


allocate(mtrx1(nmegqwan,nmegqwan))
allocate(mtrx2(nmegqwan,nmegqwan))

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
    zt1=zt1+megqwan(i,iig0q)*chi0wan(i,j)*dconjg(megqwan(j,iig0q))
    zt2=zt2+megqwan(i,iig0q)*mtrx2(i,j)*dconjg(megqwan(j,iig0q))
  enddo
enddo
f_response_(f_chi0_wann)=zt1
f_response_(f_chi_wann)=zt2
f_response_(f_epsilon_eff_wann)=1.d0/(1.d0+vhgq(iig0q,iq)*f_response_(f_chi_wann))
f_response_(f_sigma_wann)=zi*dreal(w)*(zone-f_response_(f_epsilon_eff_wann))/fourpi
f_response_(f_loss_wann)=1.d0/f_response_(f_epsilon_eff_wann)
deallocate(mtrx1)
deallocate(mtrx2)
return
end

