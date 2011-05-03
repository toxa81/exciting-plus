subroutine sic_readvwan
use modmain
use mod_sic
use mod_hdf5
use mod_linresp
use mod_wannier
implicit none
integer n,ispn,nwt,j
character*20 c1,c2
character*100 path
logical exist

inquire(file="sic.hdf5",exist=exist)
if (.not.exist) return

call hdf5_read("sic.hdf5","/","nwt",nwt)
if (nwt.ne.sic_wantran%nwt) then
  write(*,'("Error(sic_readvwan) wrong number of Wannier transitions")')
  write(*,'("  nwt read from file : ",I6)')nwt
  write(*,'("  sic_wantran%nwt : ",I6)')sic_wantran%nwt
  call pstop
endif
call hdf5_read("sic.hdf5","/","vwanme",vwanme(1),(/sic_wantran%nwt/))
call hdf5_read("sic.hdf5","/","e0",sic_wann_e0(1),(/nwantot/))
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
      call hdf5_read("sic.hdf5",path,"wanlm",s_wanlm(1,1,ispn,j),&
        (/lmmaxwan,s_nr/))
      call hdf5_read("sic.hdf5",path,"wvlm",s_wvlm(1,1,ispn,j),&
        (/lmmaxwan,s_nr/))
    enddo
  endif
enddo
tsic_wv=.true.
return
end
