subroutine sic_gndstate
use modmain
use mod_sic
implicit none
character*3 c3

sic=.true.
do isclsic=1,nsclsic
  write(c3,'(I3.3)')isclsic
  if (iproc.eq.0) write(*,'("SIC iteration : ",I4)')isclsic
  call gndstate
  call write_xml_info
  call mpi_world_barrier
  if (iproc.eq.0) then
    call system("mv INFO.OUT INFO.OUT"//c3)
    call system("mv EIGVAL.OUT EIGVAL.OUT"//c3)
    call system("mv EFERMI.OUT EFERMI.OUT"//c3)
    call system("mv TOTENERGY.OUT TOTENERGY.OUT"//c3)
    call system("mv SIC_WANN_E0.OUT SIC_WANN_E0.OUT"//c3)
    call system("mv info.xml info.xml"//c3)
  endif
  call sic_main
!  if (mpi_grid_root()) then
!    call system("cp SIC_ETOT.OUT SIC_ETOT.OUT"//c3)
!  endif
  call mpi_world_barrier
  if (iproc.eq.0) then
    call system("mv SIC.OUT SIC.OUT"//c3)  
  endif
enddo
return
end
