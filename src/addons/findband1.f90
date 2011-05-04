
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: findband
! !INTERFACE:
subroutine findband1(sol,l,k,np,nr,r,vr,de0,e,fnd,n)
! !INPUT/OUTPUT PARAMETERS:
!   sol : speed of light in atomic units (in,real)
!   l   : angular momentum quantum number (in,integer)
!   k   : quantum number k, zero if Dirac eqn. is not to be used (in,integer)
!   np  : order of predictor-corrector polynomial (in,integer)
!   nr  : number of radial mesh points (in,integer)
!   r   : radial mesh (in,real(nr))
!   vr  : potential on radial mesh (in,real(nr))
!   de0 : default energy step size (in,real)
!   e   : input energy and returned band energy (inout,real)
!   fnd : set to .true. if the band energy is found (out,logical)
! !DESCRIPTION:
!   Finds the band energies for a given radial potential and angular momentum.
!   This is done by first searching upwards in energy until the radial
!   wavefunction at the muffin-tin radius is zero. This is the energy at the top
!   of the band, denoted $E_{\rm t}$. A downward search is now performed from
!   $E_{\rm t}$ until the slope of the radial wavefunction at the muffin-tin
!   radius is zero. This energy, $E_{\rm b}$, is at the bottom of the band. The
!   band energy is taken as $(E_{\rm t}+E_{\rm b})/2$. If either $E_{\rm t}$ or
!   $E_{\rm b}$ cannot be found then the band energy is set to the input value.
!
! !REVISION HISTORY:
!   Created September 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: sol
integer, intent(in) :: l
integer, intent(in) :: k
integer, intent(in) :: np
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
real(8), intent(in) :: vr(nr)
real(8), intent(in) :: de0
real(8), intent(inout) :: e
logical, intent(out) :: fnd
integer, intent(in) :: n
! local variables
! maximum number of steps
integer maxstp
logical ttop,tbot
integer ie,nn
real(8) de,t,etop,ebot
! automatic arrays
real(8) p0(nr),p1(nr),q0(nr),q1(nr)
fnd=.false.
de=abs(de0)
etop=3.d0
maxstp=(etop+10.d0)/de
ttop=.false.
tbot=.false.
! find the top of the band
do ie=1,maxstp
  call rschroddme(sol,0,l,k,etop,np,nr,r,vr,nn,p0,p1,q0,q1)
  if (nn.eq.(n-l-1)) then
    ttop=.true.
    exit
  endif
  etop=etop-de
enddo
t=p1(nr)
ebot=etop
do ie=1,maxstp
  call rschroddme(sol,0,l,k,ebot,np,nr,r,vr,nn,p0,p1,q0,q1)
  if (p1(nr)*t.le.0) then
    tbot=.true.
    exit
  endif
  ebot=ebot-de
  t=p1(nr)
enddo
if (ttop.and.tbot) then
! set the band energy to the mid-point
  e=(etop+ebot)/2.d0
  fnd=.true.
  return
endif
! don't signal that linearization energy is not found 
!fnd=.true.
!return
end subroutine
!EOC

