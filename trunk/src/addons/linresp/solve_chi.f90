subroutine solve_chi(ngvecchi,igq0,fourpiq0,chi0m,krnl,chi_,epsilon_,lmbd)
implicit none
integer, intent(in) :: ngvecchi
integer, intent(in) :: igq0
real(8), intent(in) :: fourpiq0
complex(8), intent(in) :: chi0m(ngvecchi,ngvecchi)
complex(8), intent(in) :: krnl(ngvecchi,ngvecchi)
complex(8), intent(out) :: chi_(4)
complex(8), intent(out) :: epsilon_(4)
complex(8), intent(out) :: lmbd(ngvecchi)

complex(8), allocatable :: epsilon(:,:)
complex(8), allocatable :: mtrx1(:,:)
integer i

! Note: different epsilons and chis are introdued only to go inside the 
!   "black box" of non-linear marix equation for chi
! functions that are related to physical measurements are: 
!   1. chi_{Gq,Gq}(q-Gq,w) 
!   2. epsilon_eff(q-Gq,w)

allocate(epsilon(ngvecchi,ngvecchi))
allocate(mtrx1(ngvecchi,ngvecchi))

! save chi0
chi_(1)=chi0m(igq0,igq0)
! compute matrix 1-chi0*(v+fxc) 
epsilon=dcmplx(0.d0,0.d0)
do i=1,ngvecchi
  epsilon(i,i)=dcmplx(1.d0,0.d0)
enddo
call zgemm('N','N',ngvecchi,ngvecchi,ngvecchi,dcmplx(-1.d0,0.d0), &
  chi0m,ngvecchi,krnl,ngvecchi,dcmplx(1.d0,0.d0),epsilon,ngvecchi)
! find eigen-values of "epsilon" matrix
mtrx1=epsilon
call diagzge(ngvecchi,mtrx1,lmbd)
! save epsilon_matrix_GqGq
epsilon_(1)=epsilon(igq0,igq0)
! save epsilon_scalar_GqGq
epsilon_(2)=1.d0-chi0m(igq0,igq0)*krnl(igq0,igq0)
! invert epsilon matrix
call invzge(epsilon,ngvecchi)
! save 1/(epsilon^-1)_{GqGq}
epsilon_(3)=1.d0/epsilon(igq0,igq0)
! save chi_scalar
chi_(2)=chi0m(igq0,igq0)/epsilon_(2)
! save chi_pseudo_scalar
chi_(3)=chi0m(igq0,igq0)/epsilon_(3)
! compute chi
call zgemm('N','N',ngvecchi,ngvecchi,ngvecchi,dcmplx(1.d0,0.d0), &
  epsilon,ngvecchi,chi0m,ngvecchi,dcmplx(0.d0,0.d0),mtrx1,ngvecchi)
! save chi
chi_(4)=mtrx1(igq0,igq0)
! save epsilon_eff
epsilon_(4)=1.d0/(1.d0+fourpiq0*chi_(4))
deallocate(epsilon,mtrx1)
return
end
