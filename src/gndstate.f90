
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
integer ik,is,ia,idm,i,j,ikloc
integer n,nwork
real(8) dv,timetot
! allocatable arrays
real(8), allocatable :: v(:)
real(8), allocatable :: work(:)
real(8), allocatable :: evalfv(:,:,:)

! require forces for structural optimisation
if ((task.eq.2).or.(task.eq.3)) tforce=.true.
! initialise global variables
call timer_start(t_init,reset=.true.)
call init0
call init1
if (.not.mpi_grid_in()) return
! only root processor writes
wproc=mpi_grid_root()
! allocate arrays for eigen-values/-vectors
allocate(evalfv(nstfv,nspnfv,nkptloc))
allocate(evecfvloc(nmatmax,nstfv,nspnfv,nkptloc))
allocate(evecsvloc(nstsv,nstsv,nkptloc))
! no band disentanglement in ground state calculation
ldisentangle=.false. 
! initialise OEP variables if required
if (xctype.lt.0) call init2
if (wproc) then
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
  write(60,'("Total number of processors : ",I8)')nproc
  write(60,'("Grid dimensions : ",10I8)')mpi_grid_size  
  call writeinfo(60)
  write(60,*)
end if
! initialise or read the charge density and potentials from file
iscl=0
if ((task.eq.1).or.(task.eq.3)) then
  call readstate
  if (wproc) write(60,'("Potential read in from STATE.OUT")')
else if (task.eq.200) then
  call phveff
  if (wproc) write(60,'("Supercell potential constructed from STATE.OUT")')
else
  call rhoinit
  call poteff
  call genveffig
  if (wproc) write(60,'("Density and potential initialised from &
    &atomic data")')
end if
if (wproc) call flushifc(60)
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
if (wproc.and.((task.eq.0).or.(task.eq.2))) call delevec
call timer_stop(t_init)
! begin the self-consistent loop
if (wproc) then
  write(60,*)
  write(60,'("+------------------------------+")')
  write(60,'("| Self-consistent loop started |")')
  write(60,'("+------------------------------+")')
end if
do iscl=1,maxscl
  call timer_start(t_iter_tot,reset=.true.)
! reset all timers
  call timer_reset(t_apw_rad)
  call timer_reset(t_seceqn)
  call timer_reset(t_seceqnfv)
  call timer_reset(t_seceqnfv_setup)
  call timer_reset(t_seceqnfv_setup_h)
  call timer_reset(t_seceqnfv_setup_h_mt)
  call timer_reset(t_seceqnfv_setup_h_it)
  call timer_reset(t_seceqnfv_setup_o)
  call timer_reset(t_seceqnfv_setup_o_mt)
  call timer_reset(t_seceqnfv_setup_o_it)
  call timer_reset(t_seceqnfv_diag)
  call timer_reset(t_seceqnsv)
  call timer_reset(t_svhmlt_setup)
  call timer_reset(t_svhmlt_diag)
  call timer_reset(t_rho_mag_sum)
  call timer_reset(t_rho_mag_sym)
  call timer_reset(t_rho_mag_tot)
  call timer_reset(t_pot) 
  call timer_reset(t_dmat)   
  if (wproc) then
    write(60,*)
    write(60,'("+-------------------------+")')
    write(60,'("| Iteration number : ",I4," |")') iscl
    write(60,'("+-------------------------+")')
    call flushifc(60)
  end if
  if (iscl.ge.maxscl) then
    if (wproc) then
      write(60,*)
      write(60,'("Reached self-consistent loops maximum")')
    end if
    tlast=.true.
  end if
  if (wproc) call flushifc(60)
  call timer_start(t_apw_rad)
! generate the core wavefunctions and densities
  call gencore
! find the new linearisation energies
  call linengy
! write out the linearisation energies
  if (wproc) call writelinen
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
  call timer_stop(t_apw_rad)
! begin parallel loop over k-points
  evalsv=0.d0
  call timer_start(t_seceqn)
  do ikloc=1,nkptloc
! solve the first- and second-variational secular equations
    call seceqn(ikloc,evalfv(1,1,ikloc),evecfvloc(1,1,1,ikloc),evecsvloc(1,1,ikloc))
  end do
  call timer_stop(t_seceqn)
  call mpi_grid_reduce(evalsv(1,1),nstsv*nkpt,dims=(/dim_k/),&
    side=.true.,all=.true.)
  if (wproc) then
! find the occupation numbers and Fermi energy
    call occupy
! write out the eigenvalues and occupation numbers
    call writeeval
! write the Fermi energy to file
    call writefermi
  endif
  call mpi_grid_bcast(occsv(1,1),nstsv*nkpt,dims=(/dim_k,2/))
  if (wannier) call wann_ene_occ
! set the charge density and magnetisation to zero
  call timer_start(t_rho_mag_tot)
  call timer_start(t_rho_mag_sum)
  rhomt(:,:,:)=0.d0
  rhoir(:)=0.d0
  if (spinpol) then
    magmt(:,:,:,:)=0.d0
    magir(:,:)=0.d0
  end if
  do ikloc=1,nkptloc
! add to the density and magnetisation
    call rhomagk(ikloc,evecfvloc(1,1,1,ikloc),evecsvloc(1,1,ikloc))
  end do
  call mpi_grid_reduce(rhomt(1,1,1),lmmaxvr*nrmtmax*natmtot,dims=(/dim_k,2/),&
    all=.true.)
  call mpi_grid_reduce(rhoir(1),ngrtot,dims=(/dim_k,2/),all=.true.)
  if (spinpol) then
    call mpi_grid_reduce(magmt(1,1,1,1),lmmaxvr*nrmtmax*natmtot*ndmag,&
      dims=(/dim_k,2/),all=.true.)
    call mpi_grid_reduce(magir(1,1),ngrtot*ndmag,dims=(/dim_k,2/),&
      all=.true.)
  endif
  call rhomagsh
  call timer_stop(t_rho_mag_sum)
  call timer_start(t_rho_mag_sym)
! symmetrise the density
  call symrf(lradstp,rhomt,rhoir)
! symmetrise the magnetisation
  if (spinpol) call symrvf(lradstp,magmt,magir)
  call timer_stop(t_rho_mag_sym)
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
  call timer_stop(t_rho_mag_tot)
! LDA+U
  if (ldapu.ne.0) then
    call timer_start(t_dmat)
! generate the LDA+U density matrix
    !call gendmatlu
    call gendmatrsh
! generate the LDA+U potential matrix
    call genvmatlu
! write the LDA+U matrices to file
    if (wproc) call writeldapu
    call timer_stop(t_dmat)
  end if
! compute the effective potential
  call timer_start(t_pot)
  call poteff
  call timer_stop(t_pot)
! pack interstitial and muffin-tin effective potential and field into one array
  call packeff(.true.,n,v)
! mix in the old potential and field with the new
  call mixerifc(mixtype,n,v,dv,nwork,work)
! unpack potential and field
  call packeff(.false.,n,v)
! add the fixed spin moment effect field
  if (fixspin.ne.0) call fsmfield
! Fourier transform effective potential to G-space
  call timer_start(t_pot)
  call genveffig
  call timer_stop(t_pot)
! reduce the external magnetic fields if required
  if (reducebf.lt.1.d0) then
    bfieldc(:)=bfieldc(:)*reducebf
    bfcmt(:,:,:)=bfcmt(:,:,:)*reducebf
  end if
! compute the energy components
  call energy
  if (wproc) then
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
  end if !wproc
  call mpi_grid_bcast(tlast)
! exit self-consistent loop if last iteration is complete
  if (tlast) goto 20
! check for convergence
  if (wproc) then
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
    call timer_stop(t_iter_tot)
! output the current total CPU time
    timetot=timeinit+timemat+timefv+timesv+timerho+timepot+timefor
    write(60,*)
    write(60,'("Iteration time (seconds)                    : ",F12.2)')&
      timer_get_value(t_iter_tot)
    write(60,'("  Radial APW setup                          : ",F12.2)')&
      timer_get_value(t_apw_rad)
    write(60,'("  Total for secular equation                : ",F12.2)')&
      timer_get_value(t_seceqn)
    write(60,'("    firt-variational                        : ",F12.2)')&
      timer_get_value(t_seceqnfv)
    write(60,'("      setup                                 : ",F12.2)')&
      timer_get_value(t_seceqnfv_setup)
    write(60,'("        setup H (total, MT, IT)             : ",3F12.2)')&
      timer_get_value(t_seceqnfv_setup_h),&
      timer_get_value(t_seceqnfv_setup_h_mt),&
      timer_get_value(t_seceqnfv_setup_h_it)
    write(60,'("        setup O (total, MT, IT)             : ",3F12.2)')&
      timer_get_value(t_seceqnfv_setup_o),&
      timer_get_value(t_seceqnfv_setup_o_mt),&
      timer_get_value(t_seceqnfv_setup_o_it)
    write(60,'("      diagonalization                       : ",F12.2)')&
      timer_get_value(t_seceqnfv_diag)
    write(60,'("    second-variational                      : ",F12.2)')&
      timer_get_value(t_seceqnsv)
    write(60,'("      setup                                 : ",F12.2)')&
      timer_get_value(t_svhmlt_setup)
    write(60,'("      diagonalization                       : ",F12.2)')&
      timer_get_value(t_svhmlt_diag)
    write(60,'("  Total for charge and magnetization        : ",F12.2)')&
      timer_get_value(t_rho_mag_tot)
    write(60,'("    k-point summation                       : ",F12.2)')&
      timer_get_value(t_rho_mag_sum)
    write(60,'("    symmetrization                          : ",F12.2)')&
      timer_get_value(t_rho_mag_sym)
    write(60,'("  Total for potential                       : ",F12.2)')&
      timer_get_value(t_pot)
    write(60,'("  Density matrix setup                      : ",F12.2)')&
      timer_get_value(t_dmat)
  end if !wproc
! end the self-consistent loop

end do !iscl
20 continue
if (wproc) then
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
end if !wproc
! write eigenvalues/vectors and occupancies to file
if (mpi_grid_side(dims=(/dim_k/))) then
  do i=0,mpi_grid_size(dim_k)-1
    if (mpi_grid_x(dim_k).eq.i) then
      do ikloc=1,nkptloc
        ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
        call putevalfv(ik,evalfv(1,1,ikloc))
        call putevalsv(ik,evalsv(1,ik))
        call putevecfv(ik,evecfvloc(1,1,1,ikloc))
        call putevecsv(ik,evecsvloc(1,1,ikloc))
        call putoccsv(ik,occsv(1,ik))
!        if (wannier) call putwann(ikloc)
      end do
    end if
    call mpi_grid_barrier(dims=(/dim_k/))
  end do
endif
call mpi_grid_bcast(tstop)
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
if (wproc) then
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
