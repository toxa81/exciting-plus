
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: occupy
! !INTERFACE:
subroutine occupy2(nkpt_,wkpt_,evalsv_,occsv_)
! !USES:
use modmain
! !DESCRIPTION:
!   Finds the Fermi energy and sets the occupation numbers for the
!   second-variational states using the routine {\tt fermi}.
!
! !REVISION HISTORY:
!   Created February 2004 (JKD)
!EOP
!BOC
implicit none
integer, intent(in) :: nkpt_
real(8), intent(in) :: wkpt_(nkpt_)
real(8), intent(in) :: evalsv_(nstsv,nkpt_)
real(8), intent(out) :: occsv_(nstsv,nkpt_)
! local variables
integer, parameter :: maxit=1000
integer ik,ist,it,i,j,zval,ik0,ik1
real(8) e0,e1,chg,x,t1
real(8) etmp(nstsv,2)
! external functions
real(8) sdelta,stheta
external sdelta,stheta
! find minimum and maximum eigenvalues
e0=evalsv_(1,1)
e1=e0
do ik=1,nkpt_
  do ist=1,nstsv
    e0=min(e0,evalsv_(ist,ik))
    e1=max(e1,evalsv_(ist,ik))
  end do
end do
if (e0.lt.evalmin) then
  write(*,*)
  write(*,'("Warning(occupy2): valence eigenvalues below evalmin for s.c. &
   &loop ",I5)') iscl
end if
t1=1.d0/swidth
! determine the Fermi energy using the bisection method
do it=1,maxit
  efermi=0.5d0*(e0+e1)
  chg=0.d0
  do ik=1,nkpt_
    do ist=1,nstsv
      if (evalsv_(ist,ik).gt.evalmin) then
        x=(efermi-evalsv_(ist,ik))*t1
        occsv_(ist,ik)=occmax*stheta(stype,x)
        chg=chg+wkpt_(ik)*occsv_(ist,ik)
      else
        occsv_(ist,ik)=0.d0
      end if
    end do
  end do
  if (chg.lt.chgval) then
    e0=efermi
  else
    e1=efermi
  end if
  if ((e1-e0).lt.epsocc) goto 10
end do
write(*,*)
write(*,'("Error(occupy): could not find Fermi energy")')
write(*,*)
call pstop
10 continue
return
end subroutine
!EOC
