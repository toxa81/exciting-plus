subroutine exciting_init
use modmain
implicit none

integer ik

call readinput
call init0
call init1

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

call geturf

if (wannier) then
  do ik=1,nkpt
    call getwann(ik)
  enddo
endif

if (allocated(veffir_zfft)) deallocate(veffir_zfft)
allocate(veffir_zfft(ngrtot))
veffir_zfft(:)=veffir(:)
call zfftifc(3,ngrid,-1,veffir_zfft)

  
return
end
