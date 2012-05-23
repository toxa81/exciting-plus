subroutine solve_chi(iq,w,chi0m,krnl,f_response_)
use modmain
use mod_addons_q
use mod_expigqr
use mod_linresp
implicit none
integer, intent(in) :: iq
complex(8), intent(in) :: w
complex(8), intent(in) :: chi0m(ngq(iq),ngq(iq))
complex(8), intent(inout) :: krnl(ngq(iq),ngq(iq))
complex(8), intent(out) :: f_response_(nf_response)
! local variables
complex(8), allocatable :: epsilon(:,:)
complex(8), allocatable :: mtrx1(:,:)
complex(8), allocatable :: zm1(:,:),zm2(:,:)
integer i,ig

! Note: different epsilons and chis are introduced only to go inside the 
!   "black box" of non-linear marix equation for chi
! functions that are related to physical measurements are: 
!   1. chi_{Gq,Gq}(q-Gq,w)  -> S(q,w)
!   2. epsilon_eff(q-Gq,w)
!   3. sigma
!   4. loss

! construct full kernel
if (lrtype.eq.0) then
  do ig=1,ngq(iq)
    krnl(ig,ig)=krnl(ig,ig)+vhgq(ig,iq)
  enddo
endif
! for magnetic response
!if (lrtype.eq.1) then
!  call genixc(ixcft)
! contruct Ixc_{G,G'}=Ixc(G-G')
!  do i=1,ngq(iq)
!    do j=1,ngq(iq)
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

allocate(epsilon(ngq(iq),ngq(iq)))
allocate(mtrx1(ngq(iq),ngq(iq)))
allocate(zm1(ngq(iq),ngq(iq)))
allocate(zm2(ngq(iq),ngq(iq)))

! save chi0_GqGq
f_response_(f_chi0)=chi0m(iig0q,iig0q)
! compute matrix 1-chi0*(v+fxc) 
epsilon=zzero
do i=1,ngq(iq)
  epsilon(i,i)=zone
enddo
call zgemm('N','N',ngq(iq),ngq(iq),ngq(iq),dcmplx(-1.d0,0.d0), &
  &chi0m,ngq(iq),krnl,ngq(iq),zone,epsilon,ngq(iq))
! save epsilon_matrix_GqGq
f_response_(f_epsilon_matrix_GqGq)=epsilon(iig0q,iig0q)
! save epsilon_scalar_GqGq
f_response_(f_epsilon_scalar_GqGq)=1.d0-chi0m(iig0q,iig0q)*krnl(iig0q,iig0q)
! invert epsilon matrix
call invzge(epsilon,ngq(iq))
! save 1/(epsilon^-1)_{GqGq}
f_response_(f_inv_epsilon_inv_GqGq)=1.d0/epsilon(iig0q,iig0q)
! save (epsilon^-1)_{GqGq}
f_response_(f_epsilon_inv_GqGq)=epsilon(iig0q,iig0q)-zone
! save chi_scalar
f_response_(f_chi_scalar)=chi0m(iig0q,iig0q)/f_response_(f_epsilon_scalar_GqGq)
! save chi_pseudo_scalar
f_response_(f_chi_pseudo_scalar)=chi0m(iig0q,iig0q)/f_response_(f_epsilon_matrix_GqGq)
! compute chi=epsilon^-1 * chi0
call zgemm('N','N',ngq(iq),ngq(iq),ngq(iq),dcmplx(1.d0,0.d0), &
  &epsilon,ngq(iq),chi0m,ngq(iq),dcmplx(0.d0,0.d0),mtrx1,ngq(iq))
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
!  call zgemm('N','N',ngq(iq),ngq(iq),ngq(iq),dcmplx(1.d0,0.d0), &
!    mtrx1,ngq(iq),krnl,ngq(iq),dcmplx(0.d0,0.d0),zm2,ngq(iq))
!! compute krnl_scr=v*zm2
!  call zgemm('N','N',ngq(iq),ngq(iq),ngq(iq),dcmplx(1.d0,0.d0), &
!    krnl,ngq(iq),zm2,ngq(iq),dcmplx(1.d0,0.d0),krnl_scr,ngq(iq))
!! compute screened Coulomb potential using "symmetrized" dielectric function
!  do ig1=1,ngq(iq)
!    do ig2=1,ngq(iq)
!      epsilon(ig1,ig2)=-vcgq(ig1)*chi0m(ig1,ig2)*vcgq(ig2)
!    enddo
!    epsilon(ig1,ig1)=dcmplx(1.d0,0.d0)+epsilon(ig1,ig1)
!  enddo
!  call invzge(epsilon,ngq(iq))
!  do ig1=1,ngq(iq)
!    do ig2=1,ngq(iq)
!      zm1(ig1,ig2)=vcgq(ig1)*epsilon(ig1,ig2)*vcgq(ig2)
!    enddo
!  enddo
!! compute difference
!  d1=0.d0
!  do ig1=1,ngq(iq)
!    do ig2=1,ngq(iq)
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
