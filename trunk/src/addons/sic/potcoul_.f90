
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: potcoul
! !INTERFACE:
subroutine potcoul_(zrhomt,zrhoir,zvclmt,zvclir)
! !USES:
use modmain
! !DESCRIPTION:
!   Calculates the Coulomb potential of the real charge density stored in the
!   global variables {\tt rhomt} and {\tt rhoir} by solving Poisson's equation.
!   These variables are coverted to complex representations and passed to the
!   routine {\tt zpotcoul}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC
implicit none
complex(8), intent(in) :: zrhomt(lmmaxvr,nrmtmax,natmtot)
complex(8), intent(in) :: zrhoir(ngrtot)
complex(8), intent(out) :: zvclmt(lmmaxvr,nrmtmax,natmtot)
complex(8), intent(out) :: zvclir(ngrtot)
! local variables
integer is,ia,ias,ir,lmax
complex(8) zrho0
real(8) spzn1(maxspecies)

! allocatable arrays
real(8), allocatable :: jlgr(:,:,:)
allocate(jlgr(0:lmaxvr+npsden+1,ngvec,nspecies))
! compute the required spherical Bessel functions
lmax=lmaxvr+npsden+1
call genjlgpr(lmax,gc,jlgr)
! solve the complex Poisson's equation
spzn1=0.d0
call zpotcoul(nrmt,nrmtmax,spnrmax,spr,1,gc,jlgr,ylmg,sfacg,spzn1,zrhomt, &
 zrhoir,zvclmt,zvclir,zrho0)
deallocate(jlgr)
return
end subroutine
!EOC

