subroutine sic_gndstate
use modmain
use mod_sic
implicit none
character*3 c3

sic=.true.
do isclsic=1,nsclsic
  write(c3,'(I3.3)')isclsic
  if (iproc.eq.0) then
    write(*,*)
    write(*,'("SIC iteration : ",I4)')isclsic
  endif
  call gndstate
! save DFT variables
  xml_info%engytot=engytot
  xml_info%bandgap=bandgap
  xml_info%rws=rwigner
  xml_info%omega=omega
  xml_info%magmom(:)=mommt(1,:)
  xml_info%fftgriddens=dble(ngrtot)/omega
  xml_info%ngrtot=ngrtot
  call mpi_world_barrier
  call sic_main
  call write_xml_info
  call mpi_world_barrier
  if (iproc.eq.0) then
    call system("mv INFO.OUT INFO.OUT"//c3)
    call system("mv EIGVAL.OUT EIGVAL.OUT"//c3)
    call system("mv EFERMI.OUT EFERMI.OUT"//c3)
    call system("mv TOTENERGY.OUT TOTENERGY.OUT"//c3)
    call system("mv SIC.OUT SIC.OUT"//c3)  
    call system("mv info.xml info.xml"//c3)
  endif
enddo
return
end
