
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
integer ik,is,ia,idm,i,j
integer n,nwork
real(8) dv,timetot
! allocatable arrays
real(8), allocatable :: v(:)
real(8), allocatable :: work(:)
real(8), allocatable :: evalfv(:,:,:)

integer, external :: ikglob

! require forces for structural optimisation
if ((task.eq.2).or.(task.eq.3)) tforce=.true.
! initialise global variables
call init0
call init1

! allocate arrays for eigevvalues/vectors
allocate(evalfv(nstfv,nspnfv,nkptloc(iproc)))
allocate(evecfvloc(nmatmax,nstfv,nspnfv,nkptloc(iproc)))
allocate(evecsvloc(nstsv,nstsv,nkptloc(iproc)))

! initialise OEP variables if required
if (xctype.lt.0) call init2
if (iproc.eq.0) then
! write the real and reciprocal lattice vectors to file
  call writelat
! write interatomic distances to file
  call writeiad(.false.)
! write symmetry matrices to file
  call writesym
! output the k-point set to file
  call writekpts
! write lattice vectors and atomic positions to file
  call writegeom(.false.)
  call writenn
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
end if
! initialise or read the charge density and potentials from file
iscl=0
if ((task.eq.1).or.(task.eq.3)) then
  call readstate
  if (iproc.eq.0) write(60,'("Potential read in from STATE.OUT")')
else if (task.eq.200) then
  call phveff
  if (iproc.eq.0) write(60,'("Supercell potential constructed from STATE.OUT")')
else
  call rhoinit
  call poteff
  call genveffig
  if (iproc.eq.0) write(60,'("Density and potential initialised from &
    &atomic data")')
end if
if (iproc.eq.0) call flushifc(60)
!if (wannier.and.task.eq.1.and.maxscl.gt.1) then
!  do i=0,nproc-1
!    if (iproc.eq.i) then
!      do ik=1,nkptloc(iproc)
!        call getwann(ik)
!      end do
!    end if
!    call barrier(comm_world)
!  end do
!endif
! size of mixing vector
n=lmmaxvr*nrmtmax*natmtot+ngrtot
if (spinpol) n=n*(1+ndmag)
if (ldapu.ne.0) n=n+2*lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot
! allocate mixing arrays
allocate(v(n))
! determine the size of the mixer work array
nwork=-1
call mixerifc(mixtype,n,v,dv,nwork,v)
allocate(work(nwork))
! set stop flag
tstop=.false.
10 continue
! set last iteration flag
tlast=.false.
! delete any existing eigenvector files
if (iproc.eq.0.and.((task.eq.0).or.(task.eq.2))) call delevec
! begin the self-consistent loop
if (iproc.eq.0) then
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop started |")')
  write(60,'("+------------------------------+")')
end if
do iscl=1,maxscl
  if (iproc.eq.0) then
    write(60,*)
    write(60,'("+-------------------------+")')
    write(60,'("| Iteration number : ",I4," |")') iscl
    write(60,'("+-------------------------+")')
    call flushifc(60)
  end if
  if (iscl.ge.maxscl) then
    if (iproc.eq.0) then
      write(60,*)
      write(60,'("Reached self-consistent loops maximum")')
    end if
    tlast=.true.
  end if
  if (iproc.eq.0) call flushifc(60)
! generate the core wavefunctions and densities
  call gencore
! find the new linearisation energies
  call linengy
! write out the linearisation energies
  if (iproc.eq.0) call writelinen
! generate the APW radial functions
  call genapwfr
! generate the local-orbital radial functions
  call genlofr
! compute the overlap radial integrals
  call olprad
! compute the Hamiltonian radial integrals
  call hmlrad
  call geturf
  call genurfprod
! begin parallel loop over k-points
  evalsv=0.d0
  call timer_reset(t_seceqnfv_setup)
  call timer_reset(t_seceqnfv_diag)
  call timer_reset(t_seceqnsv_setup)
  call timer_reset(t_seceqnsv_diag)
  do ik=1,nkptloc(iproc)
! solve the first- and second-variational secular equations
    call seceqn(ik,evalfv(1,1,ik),evecfvloc(1,1,1,ik),evecsvloc(1,1,ik))
  end do
  call dsync(evalsv,nstsv*nkpt,.true.,.false.)
  if (iproc.eq.0) then
! find the occupation numbers and Fermi energy
    call occupy
! write out the eigenvalues and occupation numbers
    call writeeval
! write the Fermi energy to file
    call writefermi
  endif
  call dsync(occsv,nstsv*nkpt,.false.,.true.)
  if (wannier) call wann_ene_occ
! set the charge density and magnetisation to zero
  call timer_reset(t_rho)
  call timer_start(t_rho)
  rhomt(:,:,:)=0.d0
  rhoir(:)=0.d0
  if (spinpol) then
    magmt(:,:,:,:)=0.d0
    magir(:,:)=0.d0
  end if
  do ik=1,nkptloc(iproc)
! add to the density and magnetisation
    call rhovalk(ik,evecfvloc(1,1,1,ik),evecsvloc(1,1,ik))
  end do
  call dsync(rhomt,lmmaxvr*nrmtmax*natmtot,.true.,.true.)
  call dsync(rhoir,ngrtot,.true.,.true.)
  if (spinpol) then
    call dsync(magmt,lmmaxvr*nrmtmax*natmtot*ndmag,.true.,.true.)
    call dsync(magir,ngrtot*ndmag,.true.,.true.)
  endif
  call timer_stop(t_rho)
! symmetrise the density
  call symrf(lradstp,rhomt,rhoir)
! symmetrise the magnetisation
  if (spinpol) call symrvf(lradstp,magmt,magir)
! convert the density from a coarse to a fine radial mesh
  call rfmtctof(rhomt)
! convert the magnetisation from a coarse to a fine radial mesh
  do idm=1,ndmag
    call rfmtctof(magmt(:,:,:,idm))
  end do
! add the core density to the total density
  call addrhocr
! calculate the charges
  call charge
! calculate the moments
  if (spinpol) call moment
! normalise the density
  call rhonorm
! LDA+U
  if (ldapu.ne.0) then
! generate the LDA+U density matrix
    call gendmatlu
! generate the LDA+U potential matrix
    call genvmatlu
! write the LDA+U matrices to file
    if (iproc.eq.0) call writeldapu
    call gendmatrsh
  end if
! compute the effective potential
  call poteff
! pack interstitial and muffin-tin effective potential and field into one array
  call packeff(.true.,n,v)
! mix in the old potential and field with the new
  call mixerifc(mixtype,n,v,dv,nwork,work)
! unpack potential and field
  call packeff(.false.,n,v)
! add the fixed spin moment effect field
  if (fixspin.ne.0) call fsmfield
! Fourier transform effective potential to G-space
  call genveffig
! reduce the external magnetic fields if required
  if (reducebf.lt.1.d0) then
    bfieldc(:)=bfieldc(:)*reducebf
    bfcmt(:,:,:)=bfcmt(:,:,:)*reducebf
  end if
! compute the energy components
  call energy
  if (iproc.eq.0) then
! output energy components
    call writeengy(60)
    write(60,*)
    write(60,'("Density of states at Fermi energy : ",G18.10)') fermidos
    write(60,'(" (states/Hartree/unit cell)")')
! write total energy to TOTENERGY.OUT and flush
    write(61,'(G22.12)') engytot
    call flushifc(61)
! write DOS at Fermi energy to FERMIDOS.OUT and flush
    write(62,'(G18.10)') fermidos
    call flushifc(62)
! output charges and moments
    call writechg(60)
! write total moment to MOMENT.OUT and flush
    if (spinpol) then
      write(63,'(3G18.10)') momtot(1:ndmag)
      call flushifc(63)
    end if
! output effective fields for fixed spin moment calculations
    if (fixspin.ne.0) call writefsm(60)
! check for WRITE file
    inquire(file='WRITE',exist=exist)
    if (exist) then
      write(60,*)
      write(60,'("WRITE file exists - writing STATE.OUT")')
      call writestate
      open(50,file='WRITE')
      close(50,status='DELETE')
    end if
! write STATE.OUT file if required
    if (nwrite.ge.1) then
      if (mod(iscl,nwrite).eq.0) then
        call writestate
        write(60,*)
        write(60,'("Wrote STATE.OUT")')
      end if
    end if
  end if !iproc.eq.0
  call lsync(tlast,1,.true.)
! exit self-consistent loop if last iteration is complete
  if (tlast) goto 20
! check for convergence
  if (iproc.eq.0) then
    if (iscl.ge.2) then
      write(60,*)
      write(60,'("RMS change in effective potential (target) : ",G18.10,&
       &" (",G18.10,")")') dv,epspot
      if (dv.lt.epspot) then
        write(60,*)
        write(60,'("Potential convergence target achieved")')
        tlast=.true.
      end if
      write(65,'(G18.10)') dv
      call flushifc(65)
    end if
    if (xctype.lt.0) then
      write(60,*)
      write(60,'("Magnitude of OEP residual : ",G18.10)') resoep
    end if
! check for STOP file
    inquire(file='STOP',exist=exist)
    if (exist) then
      write(60,*)
      write(60,'("STOP file exists - stopping self-consistent loop")')
      tstop=.true.
      tlast=.true.
      open(50,file='STOP')
      close(50,status='DELETE')
    end if
! output the current total CPU time
    timetot=timeinit+timemat+timefv+timesv+timerho+timepot+timefor
    write(60,*)
    write(60,'("Time (CPU seconds) : ",F12.2)') timetot
    write(60,'("  first variational matrix setup            : ",F12.2)')timer(t_seceqnfv_setup,2)
    write(60,'("  first variational matrix diagonalization  : ",F12.2)')timer(t_seceqnfv_diag,2)
    write(60,'("  second variational matrix setup           : ",F12.2)')timer(t_seceqnsv_setup,2)
    write(60,'("  second variational matrix diagonalization : ",F12.2)')timer(t_seceqnsv_diag,2)
    write(60,'("  charge and magnetisation density setup    : ",F12.2)')timer(t_rho,2)
  end if !iproc.eq.0
! end the self-consistent loop
end do !iscl
20 continue
if (iproc.eq.0) then
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop stopped |")')
  write(60,'("+------------------------------+")')
! write density and potentials to file only if maxscl > 1
  if (maxscl.gt.1) then
    call writestate
    write(60,*)
    write(60,'("Wrote STATE.OUT")')
  end if
end if !iproc.eq.0
! write eigenvalues/vectors and occupancies to file
do i=0,nproc-1
  if (iproc.eq.i) then
    do ik=1,nkptloc(iproc)
      call putevalfv(ikglob(ik),evalfv(1,1,ik))
      call putevalsv(ikglob(ik),evalsv(1,ikglob(ik)))
      call putevecfv(ikglob(ik),evecfvloc(1,1,1,ik))
      call putevecsv(ikglob(ik),evecsvloc(1,1,ik))
      call putoccsv(ikglob(ik),occsv(1,ikglob(ik)))
!      if (wannier) call putwfc(ikglob(ik),wann_c(1,1,1,ik))
      if (wannier) call putwann(ik)
    end do
  end if
  call barrier(comm_world)
end do
if (wannier) call zsync(wann_h,nwann*nwann*nkpt,.true.,.false.)
if (wannier.and.iproc.eq.0) then
  do i=1,nwann
    wann_h(i,i,:)=wann_h(i,i,:)-efermi
  enddo
  wann_h=wann_h*ha2ev
  open(200,file='hamilt',form='formatted',status='replace')
  write(200,*)nkpt,nwann
  do ik=1,nkpt
    write(200,*)1.d0 !wtkp(ikp)
    do i=1,nwann
      do j=1,nwann
        write(200,*)dreal(wann_h(i,j,ik)),dimag(wann_h(i,j,ik))
      enddo
    enddo
  enddo	
  close(200)
endif

call lsync(tstop,1,.true.)
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
! add blank line to TOTENERGY.OUT, FERMIDOS.OUT, MOMENT.OUT and RMSDVEFF.OUT
  write(61,*)
  write(62,*)
  if (spinpol) write (63,*)
  write(65,*)
! begin new self-consistent loop with updated positions
  goto 10
end if
30 continue
! output timing information
if (iproc.eq.0) then
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
endif
deallocate(evalfv,evecfvloc,evecsvloc)
deallocate(v,work)
return
end subroutine
!EOC
