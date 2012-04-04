subroutine sic_read_data(treadk)
use modmain
use mod_sic
use mod_hdf5
use mod_linresp
use mod_wannier
implicit none
logical, intent(in) :: treadk
!
integer n,nwt,j,ik,ikloc,nkpt_
character*100 path
logical file_exists,read_vme
!
inquire(file="sic.hdf5",exist=file_exists)
tsic_wv=.false.
if (.not.file_exists) return
call hdf5_read("sic.hdf5","/","nwt",nwt)
read_vme=.true.
if (nwt.ne.sic_wantran%nwt) then
  write(*,'("Warning(sic_read_data): wrong number of Wannier transitions")')
  write(*,'("  sic.hdf5 : ",I6)')nwt
  write(*,'("  current  : ",I6)')sic_wantran%nwt
  read_vme=.false.
endif
if (read_vme) call hdf5_read("sic.hdf5","/","vme",sic_vme(1),(/sic_wantran%nwt/))
call hdf5_read("sic.hdf5","/","sic_energy_tot",sic_energy_tot)
call hdf5_read("sic.hdf5","/","sic_energy_pot",sic_energy_pot)
call hdf5_read("sic.hdf5","/","sic_energy_kin",sic_energy_kin)
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  write(path,'("/wannier_functions/",I4.4)')n
  call hdf5_read("sic.hdf5",path,"wlm",s_wlm(1,1,1,j),&
    &(/lmmaxwan,s_nr,nspinor/))
  call hdf5_read("sic.hdf5",path,"wvlm",s_wvlm(1,1,1,j),&
    &(/lmmaxwan,s_nr,nspinor/))
enddo
if (treadk) then
  call hdf5_read("sic.hdf5","/","nkpt",nkpt_)
  if (nkpt.ne.nkpt_) then
    write(*,'("Warning(sic_read_data): wrong nkpt")')
    write(*,'("  sic.hdf5 : ",I6)')nkpt_
    write(*,'("  current  : ",I6)')nkpt
    goto 10
  endif
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    do ikloc=1,nkptloc
      ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
      write(path,'("/wannier_functions/",I4.4,"/kpoints/",I4.4)')n,ik
      call hdf5_read("sic.hdf5",path,"wkmt",s_wkmt(1,1,1,1,j,ikloc),&
        &(/nrmtmax,lmmaxapw,natmtot,nspinor/))
      call hdf5_read("sic.hdf5",path,"wkit",s_wkit(1,1,j,ikloc),&
        &(/ngkmax,nspinor/))
    enddo
  enddo 
endif
10 continue
tsic_wv=.true.
return
end
