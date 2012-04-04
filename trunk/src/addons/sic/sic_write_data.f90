subroutine sic_write_data
use modmain
use mod_sic
use mod_hdf5
use mod_wannier
implicit none
integer n,j
character*100 path
!
if (wproc) then
  call hdf5_write("sic.hdf5","/","vme",sic_vme(1),(/sic_wantran%nwt/))
  call hdf5_write("sic.hdf5","/","sic_energy_tot",sic_energy_tot)
  call hdf5_write("sic.hdf5","/","sic_energy_pot",sic_energy_pot)
  call hdf5_write("sic.hdf5","/","sic_energy_kin",sic_energy_kin)  
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    write(path,'("/wannier_functions/",I4.4)')n
    call hdf5_write("sic.hdf5",path,"wlm",s_wlm(1,1,1,j),&
      &(/lmmaxwan,s_nr,nspinor/))
    call hdf5_write("sic.hdf5",path,"wvlm",s_wvlm(1,1,1,j),&
      &(/lmmaxwan,s_nr,nspinor/))
  enddo
endif
return
end
