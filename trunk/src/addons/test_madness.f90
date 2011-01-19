#ifdef _MAD_
subroutine test_madness
use modmain
use mod_nrkp
use mod_wannier
use mod_nadness
implicit none
integer i
real(8) d1,vrc(3),val(2)

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



call timer_start(40,reset=.true.)
vrc=0.d0
do i=1,10000
  vrc(1)=1.2d0
!  !call elk_wan_val(1,1,8.d0,vrc,val)
  call elk_wan_rho(1,8.d0,vrc,val)
!  write(100,*)vrc(1),val(1) !abs(val(1))**2
enddo
call timer_stop(40)
write(*,*)timer_get_value(40)
!allocate(wann_unkmt1(lmmaxvr,nufrmax,natmtot,nspinor,nwantot))
!allocate(wann_unkit1(ngkmax,nspinor,nwantot))

!call madness_init_box
!call mpi_world_barrier
!call madness_resolve_wannier(nwantot)

!call madness_genpot
!call madness_getpot(d1)
!write(*,*)d1
call bstop
return
end
#endif
