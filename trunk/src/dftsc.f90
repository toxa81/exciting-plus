subroutine dftsc(n,nu,mu,beta,f)
use        modmain
#ifdef _MPI_
use mpi
#endif
implicit   none

integer              :: n
real(8)              :: nu(n)
real(8)              :: mu(n)
real(8)              :: beta(n)
real(8)              :: f(n)

! local variables
logical exist
integer ik,is,ia,idm,ikloc,ierr,i
real(8) dv,timetot
! allocatable arrays
real(8), allocatable :: evalfv(:,:,:)

allocate(evalfv(nstfv,nspnfv,nkptloc(iproc)))
allocate(evecfvloc(nmatmax,nstfv,nspnfv,nkptloc(iproc)))
allocate(evecsvloc(nstsv,nstsv,nkptloc(iproc)))

! begin the self-consistent loop
if (iproc.eq.0) then
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop started |")')
  write(60,'("+------------------------------+")')
endif
do iscl=1,maxscl
  if (iproc.eq.0) then
    write(60,*)
    write(60,'("+-------------------------+")')
    write(60,'("| Iteration number : ",I4," |")') iscl
    write(60,'("+-------------------------+")')
    call flushifc(60)
  endif
  if (iscl.ge.maxscl) then
    if (iproc.eq.0) then
      write(60,*)
      write(60,'("Reached self-consistent loops maximum")')
    endif
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

! begin parallel loop over k-points
  evalsv=0.d0
  spnchr=0.d0
  do ikloc=1,nkptloc(iproc)
#ifdef _MPI_
    if (iproc.eq.0) write(60,'("k-step ",I4," out of ",I4)')ikloc,nkptloc(iproc)
#endif
    ik=ikptloc(iproc,1)+ikloc-1
! solve the first- and second-variational secular equations
    call seceqn(ik,evalfv(1,1,ikloc),evecfvloc(1,1,1,ikloc),evecsvloc(1,1,ikloc))
  enddo
  call dsync(evalsv,nstsv*nkpt,.true.,.false.)
  call dsync(spnchr,nspinor*nstsv*nkpt,.true.,.false.)

  if (iproc.eq.0) then 
! find the occupation numbers and Fermi energy
    call occupy
! write out the eigenvalues and occupation numbers
    call writeeval
! write the Fermi energy to file
    call writefermi
  endif
  call dsync(occsv,nstsv*nkpt,.false.,.true.)

! set the charge density and magnetisation to zero
  rhomt(:,:,:)=0.d0
  rhoir(:)=0.d0
  if (spinpol) then
    magmt(:,:,:,:)=0.d0
    magir(:,:)=0.d0
  end if
  do ikloc=1,nkptloc(iproc)
    ik=ikptloc(iproc,1)+ikloc-1
! add to the density and magnetisation
    call rhovalk(ik,evecfvloc(1,1,1,ikloc),evecsvloc(1,1,ikloc))
  end do
  call dsync(rhomt,lmmaxvr*nrmtmax*natmtot,.true.,.true.)
  call dsync(rhoir,ngrtot,.true.,.true.)
  call dsync(magmt,lmmaxvr*nrmtmax*natmtot*ndmag,.true.,.true.)
  call dsync(magir,ngrtot*ndmag,.true.,.true.)
  
  if (iproc.eq.0) then
    do ik=1,nkpt
! write the occupancies to file
      call putoccsv(ik,occsv(1,ik))
    end do
  endif
  
! symmetrise the density
  call symrf(lradstp,rhomt,rhoir)
! symmetrise the magnetisation
  if (spinpol) call symrvf(lradstp,magmt,magir)
! convert the density from a coarse to a fine radial mesh
  call rfmtctof(rhomt)
! convert the magnetisation from a coarse to a fine radial mesh
  do idm=1,ndmag
    call rfmtctof(magmt(1,1,1,idm))
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
  end if
  
! compute the effective potential
  call poteff
! pack interstitial and muffin-tin effective potential and field into one array
  call packeff(.true.,n,nu)
! mix in the old potential and field with the new
  call mixer(.false.,beta0,betamax,n,nu,mu,beta,f,dv)
! unpack potential and field
  call packeff(.false.,n,nu)
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
    write(60,'(" (states/Hartree/spin/unit cell)")')
! write total energy to TOTENERGY.OUT and flush
    write(61,'(G18.10)') engytot
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
  endif !iproc.eq.0
! exit self-consistent loop if last iteration is complete
  call lsync(tlast,1,.true.)
  if (tlast) goto 20
  
  if (iproc.eq.0) then
! check for convergence
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
      write(60,'("Magnitude of OEP residue : ",G18.10)') resoep
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
  endif !iproc.eq.0
! end the self-consistent loop
end do
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
endif !iproc.eq.0

do i=0,nproc-1
  if (iproc.eq.i) then
    do ikloc=1,nkptloc(iproc)
      ik=ikptloc(iproc,1)+ikloc-1
! write the eigenvalues/vectors to file
      call putevalfv(ik,evalfv(1,1,ikloc))
      call putevalsv(ik,evalsv(1,ik))
      call putevecfv(ik,evecfvloc(1,1,1,ikloc))
      call putevecsv(ik,evecsvloc(1,1,ikloc))
    enddo
  endif
#ifdef _MPI_
  call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif
enddo

deallocate(evalfv,evecfvloc,evecsvloc)

call lsync(tstop,1,.true.)

return  
end
