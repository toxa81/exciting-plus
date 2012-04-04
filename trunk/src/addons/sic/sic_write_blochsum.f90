subroutine sic_write_blochsum
use modmain
use mod_sic
use mod_hdf5
implicit none
!
integer i,j,ik,ikloc,n
character*100 path
!
if (mpi_grid_side(dims=(/dim_k/))) then
  do i=0,mpi_grid_dim_size(dim_k)-1
    if (mpi_grid_dim_pos(dim_k).eq.i) then
      do j=1,sic_wantran%nwan
        n=sic_wantran%iwan(j)
        do ikloc=1,nkptloc
          ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
          write(path,'("/wannier_functions/",I4.4,"/kpoints/",I4.4)')n,ik
          call hdf5_write("sic.hdf5",path,"wkmt",s_wkmt(1,1,1,1,j,ikloc),&
            &(/nrmtmax,lmmaxapw,natmtot,nspinor/))
          call hdf5_write("sic.hdf5",path,"wkit",s_wkit(1,1,j,ikloc),&
            &(/ngkmax,nspinor/))
        enddo
      enddo 
    endif
    call mpi_grid_barrier(dims=(/dim_k/))
  enddo
endif

return
end subroutine
