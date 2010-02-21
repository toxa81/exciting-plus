module mod_hdf5
use mod_mpi_grid
contains

subroutine hdf5_initialize
#ifdef _HDF5_
use hdf5
implicit none
integer ierr
call h5open_f(ierr)
#endif
end subroutine 

subroutine hdf5_create_file(fname)
#ifdef _HDF5_
use hdf5
implicit none
character(*), intent(in) :: fname
integer ierr
integer(hid_t) h5_root_id
call h5fcreate_f(trim(fname),H5F_ACC_TRUNC_F,h5_root_id,ierr)
if (ierr.ne.0) then
  write(*,*)
  write(*,'("Error(hdf5_create_file) : h5fcreate_f returned ",I6)')ierr
  write(*,'("  fname : ",A)')trim(fname)
  call pstop
endif
call h5fclose_f(h5_root_id,ierr)
if (ierr.ne.0) then
  write(*,*)
  write(*,'("Error(hdf5_create_file) : h5fclose_f returned ",I6)')ierr
  write(*,'("  fname : ",A)')trim(fname)
  call pstop
endif
#endif
end subroutine

subroutine hdf5_create_group(fname,path,gname)
#ifdef _HDF5_
use hdf5
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: gname
integer(hid_t) h5_root_id,h5_group_id,h5_new_group_id
integer ierr

call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
if (ierr.ne.0) then
  write(*,*)
  write(*,'("Error(hdf5_create_group) : h5fopen_f returned ",I6)')ierr
  write(*,'("  fname : ",A)')trim(fname)
  write(*,'("  path : ",A)')trim(path)
  write(*,'("  gname : ",A)')trim(gname)  
  call pstop
endif
call h5gopen_f(h5_root_id,trim(path),h5_group_id,ierr)
if (ierr.ne.0) then
  write(*,*)
  write(*,'("Error(hdf5_create_group) : h5gopen_f returned ",I6)')ierr
  write(*,'("  fname : ",A)')trim(fname)
  write(*,'("  path : ",A)')trim(path)
  write(*,'("  gname : ",A)')trim(gname)  
  call pstop
endif
call h5gcreate_f(h5_group_id,trim(gname),h5_new_group_id,ierr)
if (ierr.ne.0) then
  write(*,*)
  write(*,'("Error(hdf5_create_group) : h5gcreate_f returned ",I6)')ierr
  write(*,'("  fname : ",A)')trim(fname)
  write(*,'("  path : ",A)')trim(path)
  write(*,'("  gname : ",A)')trim(gname)  
  call pstop
endif
call h5gclose_f(h5_new_group_id,ierr)
if (ierr.ne.0) then
  write(*,*)
  write(*,'("Error(hdf5_create_group) : 1-st h5gclose_f returned ",I6)')ierr
  write(*,'("  fname : ",A)')trim(fname)
  write(*,'("  path : ",A)')trim(path)
  write(*,'("  gname : ",A)')trim(gname)  
  call pstop
endif
call h5gclose_f(h5_group_id,ierr)
if (ierr.ne.0) then
  write(*,*)
  write(*,'("Error(hdf5_create_group) : 2-nd h5gclose_f returned ",I6)')ierr
  write(*,'("  fname : ",A)')trim(fname)
  write(*,'("  path : ",A)')trim(path)
  write(*,'("  gname : ",A)')trim(gname)  
  call pstop
endif
call h5fclose_f(h5_root_id,ierr)
if (ierr.ne.0) then
  write(*,*)
  write(*,'("Error(hdf5_create_group) : h5fclose_f returned ",I6)')ierr
  write(*,'("  fname : ",A)')trim(fname)
  write(*,'("  path : ",A)')trim(path)
  write(*,'("  gname : ",A)')trim(gname)  
  call pstop
endif
#endif
end subroutine



end module
