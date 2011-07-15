subroutine sic_read_data(treadk)
use modmain
use mod_sic
use mod_hdf5
use mod_linresp
use mod_wannier
implicit none
logical, intent(in) :: treadk
!
integer n,ispn,nwt,j,ik,ikloc,nkpt_
character*20 c1,c2
character*100 path
logical exist
!
inquire(file="sic.hdf5",exist=exist)
if (.not.exist) return
call hdf5_read("sic.hdf5","/","nwt",nwt)
if (nwt.ne.sic_wantran%nwt) then
  write(*,'("Error(sic_read_data): wrong number of Wannier transitions")')
  write(*,'("  sic.hdf5 : ",I6)')nwt
  write(*,'("  current  : ",I6)')sic_wantran%nwt
  call pstop
endif
call hdf5_read("sic.hdf5","/","vme",sic_vme(1),(/sic_wantran%nwt/))
call hdf5_read("sic.hdf5","/","e0",sic_wan_e0(1),(/nwantot/))
call hdf5_read("sic.hdf5","/","sic_energy_tot",sic_energy_tot)
call hdf5_read("sic.hdf5","/","sic_energy_pot",sic_energy_pot)
call hdf5_read("sic.hdf5","/","sic_energy_kin",sic_energy_kin)  
do n=1,nwantot
  j=sic_wantran%idxiwan(n)
  if (j.gt.0) then
    do ispn=1,nspinor
      write(c1,'("n",I4.4)')n
      write(c2,'("s",I4.4)')ispn
      path="/wann/"//trim(adjustl(c1))//"/"//trim(adjustl(c2))
      call hdf5_read("sic.hdf5",path,"wlm",s_wlm(1,1,ispn,j),&
        (/lmmaxwan,s_nr/))
      call hdf5_read("sic.hdf5",path,"wvlm",s_wvlm(1,1,ispn,j),&
        (/lmmaxwan,s_nr/))
    enddo
  endif
enddo
if (treadk) then
  call hdf5_read("sic.hdf5","/kpoint","nkpt",nkpt_)
  if (nkpt.ne.nkpt_) then
    write(*,'("Error(sic_read_data): wrong nkpt")')
    write(*,'("  sic.hdf5 : ",I6)')nkpt_
    write(*,'("  current  : ",I6)')nkpt
    call pstop
  endif
  do ikloc=1,nkptloc
    ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
    write(c1,'(I4.4)')ik
    do n=1,nwantot
      j=sic_wantran%idxiwan(n)
      if (j.gt.0) then
        write(c2,'("n",I4.4)')n
        path="/kpoint/"//trim(adjustl(c1))//"/"//trim(adjustl(c2))
        call hdf5_read("sic.hdf5",path,"wkmt",s_wkmt(1,1,1,1,j,ikloc),&
          (/nrmtmax,lmmaxapw,natmtot,nspinor/))
        call hdf5_read("sic.hdf5",path,"wkit",s_wkit(1,1,j,ikloc),&
          (/ngkmax,nspinor/))
      endif
    enddo !n
  enddo !ikloc
endif
tsic_wv=.true.
return
end
