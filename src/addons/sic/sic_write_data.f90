subroutine sic_write_data
use modmain
use mod_sic
use mod_hdf5
use mod_wannier
implicit none
integer n,ispn,j,i,ikloc,ik
character*12 c1,c2
character*100 path
!
if (wproc) then
  call hdf5_write("sic.hdf5","/","vme",sic_vme(1),(/sic_wantran%nwt/))
  call hdf5_write("sic.hdf5","/","sic_energy_tot",sic_energy_tot)
  call hdf5_write("sic.hdf5","/","sic_energy_pot",sic_energy_pot)
  call hdf5_write("sic.hdf5","/","sic_energy_kin",sic_energy_kin)  
  do n=1,nwantot
    j=sic_wantran%idxiwan(n)
    if (j.gt.0) then
      write(path,'("/wannier_functions/",I4.4)')n
      call hdf5_write("sic.hdf5",path,"wlm",s_wlm(1,1,1,j),&
        &(/lmmaxwan,s_nr,nspinor/))
      call hdf5_write("sic.hdf5",path,"wvlm",s_wvlm(1,1,1,j),&
        &(/lmmaxwan,s_nr,nspinor/))
    endif
  enddo
endif
if (mpi_grid_side(dims=(/dim_k/))) then
  do i=0,mpi_grid_dim_size(dim_k)-1
    if (mpi_grid_dim_pos(dim_k).eq.i) then
      do n=1,nwantot
        j=sic_wantran%idxiwan(n)
        if (j.gt.0) then
          do ikloc=1,nkptloc
            ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
            write(path,'("/wannier_functions/",I4.4,"/kpoints/",I4.4)')n,ik
            call hdf5_write("sic.hdf5",path,"wkmt",s_wkmt(1,1,1,1,j,ikloc),&
              &(/nrmtmax,lmmaxapw,natmtot,nspinor/))
            call hdf5_write("sic.hdf5",path,"wkit",s_wkit(1,1,j,ikloc),&
              &(/ngkmax,nspinor/))
          enddo
        endif
      enddo 
    endif
    call mpi_grid_barrier(dims=(/dim_k/))
  enddo
endif
return
end
