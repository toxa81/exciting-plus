#ifdef _MAD_
subroutine test_madness
use modmain
use mod_nrkp
use mod_wannier
implicit none
real(8) d1

call init0
call init1
write(*,*)"Hello from Elk from process ",iproc
! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! generate the core wavefunctions and densities
call gencore
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
call getufr
call genufrp
wproc=.false.
call genwfnr(-1,.false.)

allocate(wann_unkmt1(lmmaxvr,nufrmax,natmtot,nspinor,nwantot))
allocate(wann_unkit1(ngkmax,nspinor,nwantot))

call madness_init_box
call madness_resolve_wannier(nwantot,nkptnr,nspinor)

!call madness_genpot
!call madness_getpot(d1)
!write(*,*)d1
call bstop
return
end
#endif
