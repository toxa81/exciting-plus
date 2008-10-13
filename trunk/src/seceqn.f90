
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: seceqn
subroutine seceqn(ik,evalfv,evecfv,evecsv)
! !USES:
use modmain
use modwann
! !INPUT/OUTPUT PARAMETERS:
!   ik     : k-point number (in,integer)
!   evalfv : first-variational eigenvalues (out,real(nstfv))
!   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
!   evecsv : second-variational eigenvectors (out,complex(nstsv,nstsv))
! !DESCRIPTION:
!   Solves the first- and second-variational secular equations. See routines
!   {\tt match}, {\tt seceqnfv}, {\tt seceqnss}, {\tt seceqnsv} and
!   {\tt spinchar}.
!
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(out) :: evalfv(nstfv,nspnfv)
complex(8), intent(out) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
integer ispn
logical l1
integer, external :: ikglob
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
! loop over first-variational spins (nspnfv=2 for spin-spirals only)
do ispn=1,nspnfv
! find the matching coefficients
  call match(ngk(ikglob(ik),ispn),gkc(1,ik,ispn),tpgkc(1,1,ik,ispn), &
   sfacgk(1,1,ik,ispn),apwalm(1,1,1,1,ispn))
! solve the first-variational secular equation
  call seceqnfv(nmat(ikglob(ik),ispn),ngk(ikglob(ik),ispn),igkig(1,ik,ispn),vgkc(1,1,ik,ispn), &
   apwalm(1,1,1,1,ispn),evalfv(1,ispn),evecfv(1,1,ispn))
end do
if (spinsprl) then
! solve the spin-spiral second-variational secular equation
  call seceqnss(ik,apwalm,evalfv,evecfv,evecsv)
else
! solve the second-variational secular equation
  if (wannier) then
    l1=wann_add_poco
    wann_add_poco=.false.
  endif
  call seceqnsv(ik,apwalm,evalfv,evecfv,evecsv)
  if (wannier) then                                                                                                      
    call genwfc(ik,lmaxvr,lmmaxvr,nrfmax,ufrprod,evecfv,evecsv)
    if (l1) then
      call genwfpoco(ik)
      wann_add_poco=.true.
      call seceqnsv(ik,apwalm,evalfv,evecfv,evecsv)
    endif
  endif                                                                                                                                         
end if
! compute the spin characters
call spinchar(ikglob(ik),evecsv)
deallocate(apwalm)
return
end subroutine
!EOC

