module mod_hdf5
use mod_mpi_grid

interface hdf5_write
  module procedure hdf5_write_z
end interface

interface hdf5_read
  module procedure hdf5_read_z
end interface

contains

subroutine hdf5_initialize
#ifdef _HDF5_
use hdf5
implicit none
integer ierr
call h5open_f(ierr)
#endif
end subroutine

subroutine hdf5_finalize
use hdf5
implicit none
integer ierr
call h5close_f(ierr)
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

subroutine hdf5_write_z(fname,path,dname,val,dims)
use hdf5
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
complex(8), intent(in) :: val
integer, optional, dimension(:), intent(in) :: dims

integer ndims,ierr
!integer(hsize_t), allocatable :: dims_(:)
integer, allocatable :: dims_(:)
!integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id

if (present(dims)) then
  ndims=size(dims)+1
  allocate(dims_(ndims))
  dims_(1)=2
  dims_(2:ndims)=dims(:)
else
  ndims=1
  allocate(dims_(ndims))
  dims_(1)=2
endif
call write_real8_array1(val,ndims,dims_,fname,path,dname)

!call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
!
!call h5screate_simple_f(ndims,dims_,dataspace_id,ierr)
!
!call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
!
!call h5dcreate_f(group_id,trim(dname),H5T_NATIVE_DOUBLE,dataspace_id, &
!  dataset_id,ierr)
!  
!!call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,val,dims_,ierr)
!
!call h5gclose_f(group_id,ierr)
!
!call h5sclose_f(dataspace_id,ierr)
!
!call h5dclose_f(dataset_id,ierr)
!
!call h5fclose_f(h5_root_id,ierr)
!
deallocate(dims_)

  



end subroutine

subroutine hdf5_read_z(fname,path,dname,val,dims)
use hdf5
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
complex(8), intent(out) :: val
integer, optional, dimension(:), intent(in) :: dims

integer ndims,ierr
!integer(hsize_t), allocatable :: dims_(:)
integer, allocatable :: dims_(:)
!integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id

if (present(dims)) then
  ndims=size(dims)+1
  allocate(dims_(ndims))
  dims_(1)=2
  dims_(2:ndims)=dims(:)
else
  ndims=1
  allocate(dims_(ndims))
  dims_(1)=2
endif
call read_real8_array1(val,ndims,dims_,fname,path,dname)

!call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
!
!call h5screate_simple_f(ndims,dims_,dataspace_id,ierr)
!
!call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
!
!call h5dcreate_f(group_id,trim(dname),H5T_NATIVE_DOUBLE,dataspace_id, &
!  dataset_id,ierr)
!  
!!call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,val,dims_,ierr)
!
!call h5gclose_f(group_id,ierr)
!
!call h5sclose_f(dataspace_id,ierr)
!
!call h5dclose_f(dataset_id,ierr)
!
!call h5fclose_f(h5_root_id,ierr)
!
deallocate(dims_)

  



end subroutine


end module


subroutine write_real8_array1(a,ndims,dims,fname,path,nm)
use hdf5
use mod_mpi_grid
implicit none
integer, intent(in) :: ndims
integer, intent(in) :: dims(ndims)
real(8), intent(in) :: a(*)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm

integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr,i
integer(hsize_t), dimension(ndims) :: h_dims
character*100 errmsg

do i=1,ndims
  h_dims(i)=dims(i)
enddo
call h5fopen_f(fname,H5F_ACC_RDWR_F,h5_root_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(write_real8_array) : h5fopen_f returned",I6)')ierr
  goto 100
endif
call h5screate_simple_f(ndims,h_dims,dataspace_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(write_real8_array) : h5screate_simple_f returned",I6)')ierr
  goto 100
endif
call h5gopen_f(h5_root_id,path,group_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(write_real8_array) : h5gopen_f returned",I6)')ierr
  goto 100
endif
call h5dcreate_f(group_id,nm,H5T_NATIVE_DOUBLE,dataspace_id, &
  dataset_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(write_real8_array) : h5dcreate_f returned",I6)')ierr
  goto 100
endif 
call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,a,h_dims,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(write_real8_array) : h5dwrite_f returned",I6)')ierr
  goto 100
endif 
call h5dclose_f(dataset_id,ierr)
call h5gclose_f(group_id,ierr)
call h5sclose_f(dataspace_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
100 continue
write(*,*)
write(*,'(A)')trim(errmsg)
write(*,'("  ndims : ",I4)')ndims
write(*,'("  dims : ",10I4)')dims
write(*,'("  fname : ",A)')trim(fname)
write(*,'("  path : ",A)')trim(path)
write(*,'("  nm : ",A)')trim(nm)
call pstop
end

subroutine read_real8_array1(a,ndims,dims,fname,path,nm)
use hdf5
use mod_mpi_grid
implicit none
integer, intent(in) :: ndims
integer, intent(in) :: dims(ndims)
real(8), intent(out) :: a(*)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm

integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr,i
integer(HSIZE_T), dimension(ndims) :: h_dims
character*100 errmsg


do i=1,ndims
  h_dims(i)=dims(i)
enddo

call h5fopen_f(fname,H5F_ACC_RDONLY_F,h5_root_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(read_real8_array) : h5fopen_f returned",I6)')ierr
  goto 100
endif
call h5gopen_f(h5_root_id,path,group_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(read_real8_array) : h5gopen_f returned",I6)')ierr
  goto 100
endif
call h5dopen_f(group_id,nm,dataset_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(read_real8_array) : h5dopen_f returned",I6)')ierr
  goto 100
endif
call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,a,h_dims,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(read_real8_array) : h5dread_f returned",I6)')ierr
  goto 100
endif
call h5gclose_f(group_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
100 continue
write(*,*)
write(*,'(A)')trim(errmsg)
write(*,'("  ndims : ",I4)')ndims
write(*,'("  dims : ",10I4)')dims
write(*,'("  fname : ",A)')trim(fname)
write(*,'("  path : ",A)')trim(path)
write(*,'("  nm : ",A)')trim(nm)
call pstop
end



