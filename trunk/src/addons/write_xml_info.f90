subroutine write_xml_info
use modmain
use mod_wannier
implicit none
integer ias,n
if (mpi_grid_root()) then
  open(300,file="info.xml",form="formatted",status="replace")
  write(300,'("<?xml version=""1.0""?>")')
  write(300,'("<info>")')
  do ias=1,natmtot
    write(300,'("  <atom id=""",I6,""">")')ias
    write(300,'("    <magmom units=""a.u."">",F18.10,"</magmom>")')mommt(1,ias)
    write(300,'("  </atom>")')
  enddo
  do n=1,nwantot
    write(300,'("  <wannier id=""",I6,""">")')n
    write(300,'("    <magmom units=""a.u."">",F18.10,"</magmom>")')wanmom(n)
    write(300,'("  </wannier>")')
  enddo
  write(300,'("  <rws units=""a.u."">",F18.10,"</rws>")')rwigner
  write(300,'("  <etot units=""Ha"">",F18.10,"</etot>")')engytot
  write(300,'("  <gap units=""eV"">",F18.10,"</gap>")')bandgap*ha2ev
  write(300,'("</info>")')
  close(300)
endif
return
end
