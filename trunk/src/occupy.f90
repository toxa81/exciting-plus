
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: occupy
! !INTERFACE:
subroutine occupy
! !USES:
use modmain
use modtest
! !DESCRIPTION:
!   Finds the Fermi energy and sets the occupation numbers for the
!   second-variational states using the routine {\tt fermi}.
!
! !REVISION HISTORY:
!   Created February 2004 (JKD)
!   Added gap estimation, November 2009 (F. Cricchio)
!   Added adaptive smearing width, April 2010 (T. Bjorkman)
!EOP
!BOC
implicit none
! local variables
integer, parameter :: maxit=1000
integer ik,ist,it,nval,i,j
real(8) e0,e1,e,chg,x,t1
real(8) etmp(nstsv,2)
! external functions
real(8), external :: sdelta
real(8), external :: stheta
! determine the smearing width automatically if required
if ((autoswidth).and.(iscl.gt.1)) call findswidth
! find minimum and maximum eigenvalues
e0=evalsv(1,1)
e1=e0
do ik=1,nkpt
  do ist=1,nstsv
    e=evalsv(ist,ik)
    if (e.lt.e0) e0=e
    if (e.gt.e1) e1=e
  end do
end do
!if (e0.lt.e0min-1.d0) then
!  write(*,*)
!  write(*,'("Warning(occupy): minimum eigenvalue less than minimum &
!   &linearisation energy : ",2G18.10)') e0,e0min
!  write(*,'(" for s.c. loop ",I5)') iscl
!end if
t1=1.d0/swidth
! determine the Fermi energy using the bisection method
do it=1,maxit
  efermi=0.5d0*(e0+e1)
  chg=0.d0
  do ik=1,nkpt
    do ist=1,nstsv
      x=(efermi-evalsv(ist,ik))*t1
      occsv(ist,ik)=occmax*stheta(stype,x)
      chg=chg+wkpt(ik)*occsv(ist,ik)
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
stop
10 continue
! find the density of states at the Fermi surface in units of
! states/Hartree/unit cell
fermidos=0.d0
do ik=1,nkpt
  do ist=1,nstsv
    x=(evalsv(ist,ik)-efermi)*t1
    fermidos=fermidos+wkpt(ik)*sdelta(stype,x)*t1
  end do
  if (occsv(nstsv,ik).gt.epsocc) then
    write(*,*)
    write(*,'("Warning(occupy): not enough empty states for k-point ",I6)') ik
    write(*,'(" and s.c. loop ",I5)') iscl
  end if
end do
fermidos=fermidos*occmax
! write Fermi density of states to test file
call writetest(500,'DOS at Fermi energy',tol=1.d-3,rv=fermidos)
! estimate the band gap (FC)
!e0=-1.d8
!e1=1.d8
!do ik=1,nkpt
!  do ist=1,nstsv
!    e=evalsv(ist,ik)
!    if (e.lt.efermi) then
!      if (e.gt.e0) e0=e
!    else
!      if (e.lt.e1) e1=e
!    end if
!  end do
!end do
!bandgap=e1-e0
! estimate the band gap
nval=nint(chgval)
bandgap=0.d0
if (spinpol.or.(.not.spinpol.and.mod(nval,2).eq.0)) then
  ist=nval
  if (.not.spinpol) ist=ist/2
  do j=1,nstsv
    etmp(j,1)=minval(evalsv(j,:))
    etmp(j,2)=maxval(evalsv(j,:))
  enddo
! sort bands using maximum band energies
  do i=1,nstsv-1
    do j=i+1,nstsv
      if (etmp(i,2).gt.etmp(j,2)) then
        t1=etmp(i,2)
        etmp(i,2)=etmp(j,2)
        etmp(j,2)=t1
        t1=etmp(i,1)
        etmp(i,1)=etmp(j,1)
        etmp(j,1)=t1
      endif
    enddo
  enddo
  e0=etmp(ist,2)
  e1=etmp(ist+1,1)
  if (e1.gt.e0) bandgap=e1-e0
endif
! write band gap to test file
call writetest(510,'estimated band gap',tol=1.d-2,rv=bandgap)
return
end subroutine
!EOC

