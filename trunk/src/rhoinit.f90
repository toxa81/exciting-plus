
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: rhoinit
! !INTERFACE:
subroutine rhoinit
! !USES:
use modmain
! !DESCRIPTION:
!   Initialises the crystal charge density. Inside the muffin-tins it is set to
!   the spherical atomic density. In the interstitial region it is taken to be
!   constant such that the total charge is correct. Requires that the atomic
!   densities have already been calculated.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer ir,is,ia,ias
real(8) t1,t2
! zero the charge density and magnetisation arrays
rhomt(:,:,:)=0.d0
rhoir(:)=0.d0
if (spinpol) then
  magmt(:,:,:,:)=0.d0
  magir(:,:)=0.d0
end if
t1=0.d0
t2=0.d0
do ias=1,natmtot
  is=ias2is(ias)
  ia=ias2ia(ias)
  do ir=1,nrmt(is)
    rhomt(1,ir,ias)=sprho(ir,is)/y00
    t1=t1+fourpi*sprho(ir,is)*mt_rw(ir,is)
  enddo
  t2=t2+(fourpi/3)*(rmt(is)**3)
enddo
do ir=1,ngrtot
  rhoir(ir)=(chgtot-t1)/(omega-t2)
enddo
! compute the total charge
call charge
! normalise the density
call rhonorm
return
end subroutine
!EOC

