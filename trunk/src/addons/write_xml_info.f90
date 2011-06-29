subroutine write_xml_info
use modmain
use mod_wannier
implicit none
integer ias,n
real(8) momwf

if (mpi_grid_root()) then
  open(300,file="info.xml",form="formatted",status="replace")
  write(300,'("<?xml version=""1.0""?>")')
  write(300,'("<info>")')
  do ias=1,natmtot
    momwf=0.d0
    do n=1,nwantot
      if (wan_info(1,n).eq.ias) momwf=momwf+wanmom(n)
    enddo
    write(300,'("  <atom id=""",I6,""">")')ias
    write(300,'("    <magmom units=""a.u."">",F18.10,"</magmom>")')xml_info%magmom(ias)
    write(300,'("    <magmom_wf units=""a.u."">",F18.10,"</magmom_wf>")')momwf
    write(300,'("  </atom>")')
  enddo
  do n=1,nwantot
    write(300,'("  <wannier id=""",I6,""">")')n
    write(300,'("    <magmom units=""a.u."">",F18.10,"</magmom>")')wanmom(n)
    write(300,'("    <spread units=""a.u.^2"">",F18.10,"</spread>")')xml_info%wan_spread(n)
    write(300,'("  </wannier>")')
  enddo
  write(300,'("  <rws units=""a.u."">",F18.10,"</rws>")')xml_info%rws
  write(300,'("  <omega units=""a.u.^3"">",F18.10,"</omega>")')xml_info%omega
  write(300,'("  <fftgriddens units=""1/a.u.^3"">",F18.10,"</fftgriddens>")')xml_info%fftgriddens
  write(300,'("  <ngrtot>",I8,"</ngrtot>")')xml_info%ngrtot
  write(300,'("  <engytot units=""Ha"">",F18.10,"</engytot>")')xml_info%engytot
  write(300,'("  <bandgap units=""eV"">",F18.10,"</bandgap>")')xml_info%bandgap*ha2ev
  write(300,'("  <wan_tot_spread units=""a.u.^2"">",F18.10,"</wan_tot_spread>")')&
    xml_info%wan_tot_spread
  write(300,'("  <sic_energy_tot units=""Ha"">",F18.10,"</sic_energy_tot>")')&
    xml_info%sic_energy_tot
  write(300,'("  <sic_energy_pot units=""Ha"">",F18.10,"</sic_energy_pot>")')&
    xml_info%sic_energy_pot
  write(300,'("  <sic_energy_kin units=""Ha"">",F18.10,"</sic_energy_kin>")')&
    xml_info%sic_energy_kin
  write(300,'("  <sic_vme_rms>",F18.10,"</sic_vme_rms>")')xml_info%sic_vme_rms
  write(300,'("  <sic_vme_err>",F18.10,"</sic_vme_err>")')xml_info%sic_vme_err  
  write(300,'("</info>")')
  close(300)
endif
return
end
