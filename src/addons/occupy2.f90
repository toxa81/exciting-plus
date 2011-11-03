
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
integer ik,ist,it
real(8) e0,e1,chg,x,t1
! external functions
real(8) sdelta,stheta
external sdelta,stheta
! determine the smearing width automatically if required
if ((autoswidth).and.(iscl.gt.1)) call findswidth
! find minimum and maximum eigenvalues
e0=minval(evalsv_(:,:))
e1=maxval(evalsv_(:,:))
t1=1.d0/swidth
! estimate band gap
call getbandgap(nkpt_,evalsv_,bandgap,efermi)
if (bandgap.gt.0.d0) then
  do ik=1,nkpt_
    do ist=1,nstsv
      if (evalsv_(ist,ik).lt.efermi) then
        occsv_(ist,ik)=occmax
      else
        occsv_(ist,ik)=0.d0
      endif
    enddo
  enddo
else
! determine the Fermi energy using the bisection method
  do it=1,maxit
    efermi=0.5d0*(e0+e1)
    chg=0.d0
    do ik=1,nkpt_
      do ist=1,nstsv
        x=(efermi-evalsv_(ist,ik))*t1
        occsv_(ist,ik)=occmax*stheta(stype,x)
        chg=chg+wkpt_(ik)*occsv_(ist,ik)
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
  write(*,'("Error(occupy2): could not find Fermi energy")')
  write(*,*)
  call pstop
endif
10 continue
return
end subroutine
!EOC
