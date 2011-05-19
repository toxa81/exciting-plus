subroutine sic_writevwan
use modmain
use mod_sic
use mod_hdf5
use mod_linresp
use mod_wannier
implicit none
integer n,ispn,j
character*12 c1,c2
character*100 path
!
if (wproc) then
  call hdf5_create_file("sic.hdf5")
  call hdf5_create_group("sic.hdf5","/","wann")
  do n=1,nwantot
    path="/wann"
    write(c1,'("n",I4.4)')n
    call hdf5_create_group("sic.hdf5",path,trim(adjustl(c1)))   
    do ispn=1,nspinor
      path="/wann/"//trim(adjustl(c1))
      write(c2,'("s",I4.4)')ispn
      call hdf5_create_group("sic.hdf5",path,trim(adjustl(c2)))   
    enddo
  enddo
  call hdf5_write("sic.hdf5","/","nwt",sic_wantran%nwt)
  call hdf5_write("sic.hdf5","/","iwt",sic_wantran%iwt(1,1),&
    (/5,sic_wantran%nwt/))
  call hdf5_write("sic.hdf5","/","vme",sic_vme(1),(/sic_wantran%nwt/))
  call hdf5_write("sic.hdf5","/","e0",sic_wan_e0(1),(/nwantot/))
  call hdf5_write("sic.hdf5","/","sic_energy_tot",sic_energy_tot)
  call hdf5_write("sic.hdf5","/","sic_energy_pot",sic_energy_pot)
  call hdf5_write("sic.hdf5","/","sic_energy_kin",sic_energy_kin)  
  do n=1,nwantot
    j=sic_wantran%idxiwan(n)
    if (j.gt.0) then
      do ispn=1,nspinor
        write(c1,'("n",I4.4)')n
        write(c2,'("s",I4.4)')ispn
        path="/wann/"//trim(adjustl(c1))//"/"//trim(adjustl(c2))       
        call hdf5_write("sic.hdf5",path,"wlm",s_wlm(1,1,ispn,j),&
          (/lmmaxwan,s_nr/))
        call hdf5_write("sic.hdf5",path,"wvlm",s_wvlm(1,1,ispn,j),&
          (/lmmaxwan,s_nr/))
      enddo
    endif
  enddo
endif
return
end
