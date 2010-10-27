subroutine sic_gndstate
use modmain
implicit none
character*3 c3

sic=.true.
do isclsic=1,nsclsic
  write(c3,'(I3.3)')isclsic
  if (mpi_grid_root()) then
    write(*,'("SIC iteration : ",I4)')isclsic
  endif
  call gndstate
  call mpi_world_barrier
  if (mpi_grid_root()) then
    call system("cp INFO.OUT INFO.OUT"//c3)
    call system("cp EIGVAL.OUT EIGVAL.OUT"//c3)
    call system("cp EFERMI.OUT EFERMI.OUT"//c3)
    call system("cp TOTENERGY.OUT TOTENERGY.OUT"//c3)    
  endif
  call sic_genvwan
  call mpi_world_barrier
  if (mpi_grid_root()) then
    call system("cp SIC.OUT SIC.OUT"//c3)  
  endif
enddo

return
end