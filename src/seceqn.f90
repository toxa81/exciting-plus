
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: seceqn
subroutine seceqn(ikloc,evalfv,evecfv,evecsv)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ik     : k-point number (in,integer)
!   evalfv : first-variational eigenvalues (out,real(nstfv))
!   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
!   evecsv : second-variational eigenvectors (out,complex(nstsv,nstsv))
! !DESCRIPTION:
!   Solves the first- and second-variational secular equations. See routines
!   {\tt match}, {\tt seceqnfv}, {\tt seceqnss} and {\tt seceqnsv}.
!
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ikloc
real(8), intent(out) :: evalfv(nstfv,nspnfv)
complex(8), intent(out) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
integer ispn,ik,ist
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
! loop over first-variational spins (nspnfv=2 for spin-spirals only)
do ispn=1,nspnfv
! find the matching coefficients
  call match(ngk(ispn,ik),gkc(:,ispn,ikloc),tpgkc(:,:,ispn,ikloc), &
   sfacgk(:,:,ispn,ikloc),apwalm(:,:,:,:,ispn))
! solve the first-variational secular equation
  if (tseqit) then
! iteratively
    call seceqnit(nmat(ispn,ik),ngk(ispn,ik),igkig(:,ispn,ikloc),vkl(:,ik), &
     vgkl(:,:,ispn,ikloc),vgkc(:,:,ispn,ikloc),apwalm(:,:,:,:,ispn), &
     evalfv(:,ispn),evecfv(:,:,ispn))
  else
! directly
    call seceqnfv(ik,nmat(ispn,ik),ngk(ispn,ik),igkig(:,ispn,ikloc), &
     vgkc(:,:,ispn,ikloc),apwalm(:,:,:,:,ispn),evalfv(:,ispn),evecfv(:,:,ispn))
  end if
end do
if (spinsprl) then
! solve the spin-spiral second-variational secular equation
  call seceqnss(ikloc,apwalm,evalfv,evecfv,evecsv)
else
! solve the second-variational secular equation
  call seceqnsv(ikloc,apwalm,evalfv,evecfv,evecsv)
end if
! apply scissor correction if required
if (scissor.ne.0.d0) then
  do ist=1,nstsv
    if (evalsv(ist,ik).gt.efermi) evalsv(ist,ik)=evalsv(ist,ik)+scissor
  end do
end if
if (wannier) then
  call genwann(ikloc,evecfv,evecsv)
  if (wann_add_poco) then
    call wann_seceqn(ikloc,evecsv)
    call genwann(ikloc,evecfv,evecsv)
  endif
endif
deallocate(apwalm)
return
end subroutine
!EOC

