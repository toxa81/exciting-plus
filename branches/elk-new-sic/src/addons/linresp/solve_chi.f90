subroutine solve_chi(iq,w,chi0m,krnl,f_response_)
use modmain
use mod_addons_q
implicit none
integer, intent(in) :: iq
complex(8), intent(in) :: w
complex(8), intent(in) :: chi0m(ngvecme,ngvecme)
complex(8), intent(inout) :: krnl(ngvecme,ngvecme)
!complex(8), intent(out) :: krnl_scr(ngvecme,ngvecme)
complex(8), intent(out) :: f_response_(nf_response)
! local variables
complex(8), allocatable :: epsilon(:,:)
complex(8), allocatable :: mtrx1(:,:)
complex(8), allocatable :: zm1(:,:),zm2(:,:)
real(8) d1
integer i,ig1,ig2,ig

! Note: different epsilons and chis are introduced only to go inside the 
!   "black box" of non-linear marix equation for chi
! functions that are related to physical measurements are: 
!   1. chi_{Gq,Gq}(q-Gq,w)  -> S(q,w)
!   2. epsilon_eff(q-Gq,w)
!   3. sigma
!   4. loss

! construct full kernel
if (lrtype.eq.0) then
  do ig=1,ngvecme
    krnl(ig,ig)=krnl(ig,ig)+vhgq(ig,iq)
  enddo
endif
! for magnetic response
!if (lrtype.eq.1) then
!  call genixc(ixcft)
! contruct Ixc_{G,G'}=Ixc(G-G')
!  do i=1,ngvecme
!    do j=1,ngvecme
!      iv(:)=-ivg(:,gvecchi1+i-1)+ivg(:,gvecchi1+j-1)
!      krnl_rpa(i,j)=ixcft(ivgig(iv(1),iv(2),iv(3)))
!    enddo
!  enddo
!endif !lrtype.eq.1
if (lrtype.eq.1) then
  write(*,*)
  write(*,'("Error(solve_chi): Ixc kernel is required for magnetic response")')
  write(*,'("  not yet implemented")')
  write(*,*)
  call pstop
endif

allocate(epsilon(ngvecme,ngvecme))
allocate(mtrx1(ngvecme,ngvecme))
allocate(zm1(ngvecme,ngvecme))
allocate(zm2(ngvecme,ngvecme))

! save chi0_GqGq
f_response_(f_chi0)=chi0m(iig0q,iig0q)
! compute matrix 1-chi0*(v+fxc) 
epsilon=zzero
do i=1,ngvecme
  epsilon(i,i)=zone
enddo
call zgemm('N','N',ngvecme,ngvecme,ngvecme,dcmplx(-1.d0,0.d0), &
  chi0m,ngvecme,krnl,ngvecme,zone,epsilon,ngvecme)
! save epsilon_matrix_GqGq
f_response_(f_epsilon_matrix_GqGq)=epsilon(iig0q,iig0q)
! save epsilon_scalar_GqGq
f_response_(f_epsilon_scalar_GqGq)=1.d0-chi0m(iig0q,iig0q)*krnl(iig0q,iig0q)
! invert epsilon matrix
call invzge(epsilon,ngvecme)
! save 1/(epsilon^-1)_{GqGq}
f_response_(f_inv_epsilon_inv_GqGq)=1.d0/epsilon(iig0q,iig0q)
! save chi_scalar
f_response_(f_chi_scalar)=chi0m(iig0q,iig0q)/f_response_(f_epsilon_scalar_GqGq)
! save chi_pseudo_scalar
f_response_(f_chi_pseudo_scalar)=chi0m(iig0q,iig0q)/f_response_(f_epsilon_matrix_GqGq)
! compute chi=epsilon^-1 * chi0
call zgemm('N','N',ngvecme,ngvecme,ngvecme,dcmplx(1.d0,0.d0), &
  epsilon,ngvecme,chi0m,ngvecme,dcmplx(0.d0,0.d0),mtrx1,ngvecme)
! save chi
f_response_(f_chi)=mtrx1(iig0q,iig0q)
! save epsilon_eff
f_response_(f_epsilon_eff)=1.d0/(1.d0+vhgq(iig0q,iq)*f_response_(f_chi))
! save epsilon_eff_scalar
f_response_(f_epsilon_eff_scalar)=1.d0/(1.d0+vhgq(iig0q,iq)*f_response_(f_chi_scalar))

f_response_(f_sigma)=zi*dreal(w)*(zone-f_response_(f_epsilon_eff))/fourpi
f_response_(f_sigma_scalar)=zi*dreal(w)*(zone-f_response_(f_epsilon_eff_scalar))/fourpi
f_response_(f_loss)=1.d0/f_response_(f_epsilon_eff)
f_response_(f_loss_scalar)=1.d0/f_response_(f_epsilon_eff_scalar)



!if (screened_w) then
!! compute screened Coulomb potential: vscr=vbare+vbare*chi*vbare
!  krnl_scr=krnl
!! compute zm2=chi*v
!  call zgemm('N','N',ngvecme,ngvecme,ngvecme,dcmplx(1.d0,0.d0), &
!    mtrx1,ngvecme,krnl,ngvecme,dcmplx(0.d0,0.d0),zm2,ngvecme)
!! compute krnl_scr=v*zm2
!  call zgemm('N','N',ngvecme,ngvecme,ngvecme,dcmplx(1.d0,0.d0), &
!    krnl,ngvecme,zm2,ngvecme,dcmplx(1.d0,0.d0),krnl_scr,ngvecme)
!! compute screened Coulomb potential using "symmetrized" dielectric function
!  do ig1=1,ngvecme
!    do ig2=1,ngvecme
!      epsilon(ig1,ig2)=-vcgq(ig1)*chi0m(ig1,ig2)*vcgq(ig2)
!    enddo
!    epsilon(ig1,ig1)=dcmplx(1.d0,0.d0)+epsilon(ig1,ig1)
!  enddo
!  call invzge(epsilon,ngvecme)
!  do ig1=1,ngvecme
!    do ig2=1,ngvecme
!      zm1(ig1,ig2)=vcgq(ig1)*epsilon(ig1,ig2)*vcgq(ig2)
!    enddo
!  enddo
!! compute difference
!  d1=0.d0
!  do ig1=1,ngvecme
!    do ig2=1,ngvecme
!      d1=d1+abs(krnl_scr(ig1,ig2)-zm1(ig1,ig2))
!    enddo
!  enddo
!  if (d1.gt.1d-6) then
!    write(*,*)
!    write(*,'("Error(solve_chi): screened kernels must be identical")')
!    write(*,'("  difference : ",G18.10)')d1
!    call pstop
!  endif
!endif
deallocate(epsilon,mtrx1,zm1,zm2)
return
end