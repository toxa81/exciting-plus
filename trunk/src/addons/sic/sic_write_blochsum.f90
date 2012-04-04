subroutine sic_write_blochsum
use modmain
use mod_sic
use mod_hdf5
implicit none
!
integer i,j,ik,ikloc,n
character*10 c1
character*100 path
!
if (mpi_grid_root()) then
  call hdf5_create_file("sic_wnk.hdf5")
  call hdf5_create_group("sic_wnk.hdf5","/","wannier_functions")
  call hdf5_write("sic_wnk.hdf5","/","nkpt",nkpt)
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    path="/wannier_functions"
    write(c1,'(I4.4)')n
    call hdf5_create_group("sic_wnk.hdf5",path,trim(adjustl(c1)))
    path=trim(path)//"/"//trim(adjustl(c1))
    call hdf5_create_group("sic_wnk.hdf5",path,"kpoints") 
    path=trim(path)//"/"//"kpoints"
    do ik=1,nkpt
      write(c1,'(I4.4)')ik
      call hdf5_create_group("sic_wnk.hdf5",path,trim(adjustl(c1)))
    enddo
  enddo
endif
call mpi_grid_barrier()
if (mpi_grid_side(dims=(/dim_k/))) then
  do i=0,mpi_grid_dim_size(dim_k)-1
    if (mpi_grid_dim_pos(dim_k).eq.i) then
      do j=1,sic_wantran%nwan
        n=sic_wantran%iwan(j)
        do ikloc=1,nkptloc
          ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
          write(path,'("/wannier_functions/",I4.4,"/kpoints/",I4.4)')n,ik
          call hdf5_write("sic_wnk.hdf5",path,"wkmt",s_wkmt(1,1,1,1,j,ikloc),&
            &(/nrmtmax,lmmaxapw,natmtot,nspinor/))
          call hdf5_write("sic_wnk.hdf5",path,"wkit",s_wkit(1,1,j,ikloc),&
            &(/ngkmax,nspinor/))
        enddo
      enddo 
    endif
    call mpi_grid_barrier(dims=(/dim_k/))
  enddo
endif

return
end subroutine
