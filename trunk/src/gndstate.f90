
! Copyright (C) 2002-2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gndstate
! !INTERFACE:
subroutine gndstate
! !USES:
use modmain
! !DESCRIPTION:
!   Computes the self-consistent Kohn-Sham ground-state. General information is
!   written to the file {\tt INFO.OUT}. First- and second-variational
!   eigenvalues, eigenvectors and occupancies are written to the unformatted
!   files {\tt EVALFV.OUT}, {\tt EVALSV.OUT}, {\tt EVECFV.OUT}, {\tt EVECSV.OUT}
!   and {\tt OCCSV.OUT}.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!EOP
!BOC
implicit none
! local variables
logical exist
integer ik,is,ia,idm,n,i
real(8) dv,timetot
! allocatable arrays
real(8), allocatable :: nu(:)
real(8), allocatable :: mu(:)
real(8), allocatable :: beta(:)
real(8), allocatable :: f(:)
character*100 arg1

! require forces for structural optimisation
if ((task.eq.2).or.(task.eq.3)) tforce=.true.
! initialise global variables
call init0
call init1
if (wannier) then
  call wann_init
endif

! initialise OEP variables if required
if (xctype.lt.0) call init2
if (iproc.eq.0) then
! write the real and reciprocal lattice vectors to file
  call writelat
! write inter-atomic distances to file
  call writeiad
! write symmetry matrices to file
  call writesym
! output the k-point set to file
  call writekpts
! write lattice vectors and atomic positions to file
  call writegeom(.false.)
! write nearest neighbours
  call writenn
! write G-shells
  call writegshells
! open INFO.OUT file
  open(60,file='INFO'//trim(filext),action='WRITE',form='FORMATTED')
! open TOTENERGY.OUT
  open(61,file='TOTENERGY'//trim(filext),action='WRITE',form='FORMATTED')
! open FERMIDOS.OUT
  open(62,file='FERMIDOS'//trim(filext),action='WRITE',form='FORMATTED')
! open MOMENT.OUT if required
  if (spinpol) open(63,file='MOMENT'//trim(filext),action='WRITE', &
    form='FORMATTED')
! open FORCEMAX.OUT if required
  if (tforce) open(64,file='FORCEMAX'//trim(filext),action='WRITE', &
    form='FORMATTED')
! open RMSDVEFF.OUT
  open(65,file='RMSDVEFF'//trim(filext),action='WRITE',form='FORMATTED')
! write out general information to INFO.OUT
  write(60,'("Running on ",I5," proc")')nproc
  call writeinfo(60)
  write(60,*)
  call getarg(1,arg1)
  if (trim(arg1).eq.'info') call pstop
endif !iproc.eq.0
! initialise or read the charge density and potentials from file
iscl=0
if ((task.eq.1).or.(task.eq.3)) then
  call readstate
  if (iproc.eq.0) write(60,'("Density and potential read in from STATE.OUT")')
else
  call rhoinit
  call poteff
  call genveffig
  if (iproc.eq.0) write(60,'("Density and potential initialised from atomic data")')
end if
if (iproc.eq.0) call flushifc(60)
! size of mixing vector
n=lmmaxvr*nrmtmax*natmtot+ngrtot
if (spinpol) n=n*(1+ndmag)
if (ldapu.ne.0) n=n+2*lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot
! allocate mixing arrays
allocate(nu(n))
allocate(mu(n))
allocate(beta(n))
allocate(f(n))

! set stop flag
tstop=.false.
10 continue
! set last iteration flag
tlast=.false.
! initialise the mixer
call packeff(.true.,n,nu)
call mixer(.true.,beta0,betamax,n,nu,mu,beta,f,dv)
call packeff(.false.,n,nu)
! delete any existing eigenvector files
if (iproc.eq.0) call delevec

! main DFT SC loop
call dftsc(n,nu,mu,beta,f)

#ifdef _MPI_
if (tforce) then
  write(*,*)'Forces not implemented in parallel mode'
  call pstop
endif
#endif

!-----------------------!
!     compute forces    !
!-----------------------!
if ((.not.tstop).and.(tforce)) then
  call force
! output forces to INFO.OUT
  call writeforce(60)
! write maximum force magnitude to FORCEMAX.OUT
  write(64,'(G18.10)') forcemax
  call flushifc(64)
end if
!---------------------------------------!
!     perform structural relaxation     !
!---------------------------------------!
if ((.not.tstop).and.((task.eq.2).or.(task.eq.3))) then
  write(60,*)
  write(60,'("Maximum force magnitude (target) : ",G18.10," (",G18.10,")")') &
   forcemax,epsforce
  call flushifc(60)
! check force convergence
  if (forcemax.le.epsforce) then
    write(60,*)
    write(60,'("Force convergence target achieved")')
    goto 30
  end if
! update the atomic positions if forces are not converged
  call updatpos
  write(60,*)
  write(60,'("+--------------------------+")')
  write(60,'("| Updated atomic positions |")')
  write(60,'("+--------------------------+")')
  do is=1,nspecies
    write(60,*)
    write(60,'("Species : ",I4," (",A,")")') is,trim(spsymb(is))
    write(60,'(" atomic positions (lattice) :")')
    do ia=1,natoms(is)
      write(60,'(I4," : ",3F14.8)') ia,atposl(:,ia,is)
    end do
  end do
! add blank line to TOTENERGY.OUT, FERMIDOS.OUT and MOMENT.OUT
  write(61,*)
  write(62,*)
  if (spinpol) write (63,*)
! begin new self-consistent loop with updated positions
  goto 10
end if
30 continue
if (iproc.eq.0) then
! output timing information
  write(60,*)
  write(60,'("Timings (CPU seconds) :")')
  write(60,'(" initialisation",T40,": ",F12.2)') timeinit
  write(60,'(" Hamiltonian and overlap matrix set up",T40,": ",F12.2)') timemat
  write(60,'(" first-variational secular equation",T40,": ",F12.2)') timefv
  if (spinpol) then
    write(60,'(" second-variational calculation",T40,": ",F12.2)') timesv
  end if
  write(60,'(" charge density calculation",T40,": ",F12.2)') timerho
  write(60,'(" potential calculation",T40,": ",F12.2)') timepot
  if (tforce) then
    write(60,'(" force calculation",T40,": ",F12.2)') timefor
  end if
  timetot=timeinit+timemat+timefv+timesv+timerho+timepot+timefor
  write(60,'(" total",T40,": ",F12.2)') timetot
  write(60,*)
  write(60,'("+----------------------------------+")')
  write(60,'("| EXCITING version ",I1.1,".",I1.1,".",I3.3," stopped |")') version
  write(60,'("+----------------------------------+")')
! close the INFO.OUT file
  close(60)
! close the TOTENERGY.OUT file
  close(61)
! close the FERMIDOS.OUT file
  close(62)
! close the MOMENT.OUT file
  if (spinpol) close(63)
! close the FORCEMAX.OUT file
  if (tforce) close(64)
! close the RMSDVEFF.OUT file
  close(65)
endif !iproc.eq.0

deallocate(nu,mu,beta,f,nkptloc,ikptloc)

return
end subroutine
!EOC
