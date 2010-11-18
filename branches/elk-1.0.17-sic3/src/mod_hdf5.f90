module mod_hdf5
#ifdef _MPI_
use mod_mpi_grid
#endif

interface hdf5_write
  module procedure hdf5_write_z,hdf5_write_i,hdf5_write_d
end interface

interface hdf5_read
  module procedure hdf5_read_z,hdf5_read_i,hdf5_read_d
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
#ifdef _HDF5_
use hdf5
implicit none
integer ierr
call h5close_f(ierr)
#endif
end subroutine

subroutine hdf5_create_file(fname)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
#ifdef _HDF5_
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

!--------------------------!
!     hdf5_create_group    !
!--------------------------!
subroutine hdf5_create_group(fname,path,gname)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: gname
#ifdef _HDF5_
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

!---------------------!
!     hdf5_write_z    !
!---------------------!
subroutine hdf5_write_z(fname,path,dname,val,dims)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
complex(8), intent(in) :: val
integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
integer ndims
integer, allocatable :: dims_(:)
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
call write_real8_array(val,ndims,dims_,fname,path,dname)
deallocate(dims_)
#endif
end subroutine

!---------------------!
!     hdf5_write_i    !
!---------------------!
subroutine hdf5_write_i(fname,path,dname,val,dims)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
integer, intent(in) :: val
integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
integer ndims
integer, allocatable :: dims_(:)
if (present(dims)) then
  ndims=size(dims)
  allocate(dims_(ndims))
  dims_(1:ndims)=dims(:)
else
  ndims=1
  allocate(dims_(ndims))
  dims_(1)=1
endif
call write_integer_array(val,ndims,dims_,fname,path,dname)
deallocate(dims_)
#endif
end subroutine

!---------------------!
!     hdf5_write_d    !
!---------------------!
subroutine hdf5_write_d(fname,path,dname,val,dims)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
real(8), intent(in) :: val
integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
integer ndims
integer, allocatable :: dims_(:)
if (present(dims)) then
  ndims=size(dims)
  allocate(dims_(ndims))
  dims_(1:ndims)=dims(:)
else
  ndims=1
  allocate(dims_(ndims))
  dims_(1)=1
endif
call write_real8_array(val,ndims,dims_,fname,path,dname)
deallocate(dims_)
#endif
end subroutine

!--------------------!
!     hdf5_read_z    !
!--------------------!
subroutine hdf5_read_z(fname,path,dname,val,dims)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
complex(8), intent(out) :: val
integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
integer ndims
integer, allocatable :: dims_(:)
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
call read_real8_array(val,ndims,dims_,fname,path,dname)
deallocate(dims_)
#endif
end subroutine

!--------------------!
!     hdf5_read_i    !
!--------------------!
subroutine hdf5_read_i(fname,path,dname,val,dims)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
integer, intent(out) :: val
integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
integer ndims
integer, allocatable :: dims_(:)
if (present(dims)) then
  ndims=size(dims)
  allocate(dims_(ndims))
  dims_(1:ndims)=dims(:)
else
  ndims=1
  allocate(dims_(ndims))
  dims_(1)=1
endif
call read_integer_array(val,ndims,dims_,fname,path,dname)
deallocate(dims_)
#endif
end subroutine

!--------------------!
!     hdf5_read_d    !
!--------------------!
subroutine hdf5_read_d(fname,path,dname,val,dims)
#ifdef _HDF5_
use hdf5
#endif
implicit none
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: dname
real(8), intent(out) :: val
integer, optional, dimension(:), intent(in) :: dims
#ifdef _HDF5_
integer ndims
integer, allocatable :: dims_(:)
if (present(dims)) then
  ndims=size(dims)
  allocate(dims_(ndims))
  dims_(1:ndims)=dims(:)
else
  ndims=1
  allocate(dims_(ndims))
  dims_(1)=1
endif
call read_real8_array(val,ndims,dims_,fname,path,dname)
deallocate(dims_)
#endif
end subroutine

end module




#ifdef _HDF5_
!---------------------------!
!     write_real8_array     !
!---------------------------!
subroutine write_real8_array(a,ndims,dims,fname,path,nm)
#ifdef _MPI_
use mod_mpi_grid
#endif
use hdf5
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
call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(write_real8_array) : h5fopen_f returned",I6)')ierr
  goto 100
endif
call h5screate_simple_f(ndims,h_dims,dataspace_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(write_real8_array) : h5screate_simple_f returned",I6)')ierr
  goto 100
endif
call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(write_real8_array) : h5gopen_f returned",I6)')ierr
  goto 100
endif
call h5dcreate_f(group_id,trim(nm),H5T_NATIVE_DOUBLE,dataspace_id, &
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

!--------------------------!
!     read_real8_array     !
!--------------------------!
subroutine read_real8_array(a,ndims,dims,fname,path,nm)
#ifdef _MPI_
use mod_mpi_grid
#endif
use hdf5
implicit none
integer, intent(in) :: ndims
integer, intent(in) :: dims(ndims)
real(8), intent(out) :: a(*)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm

integer(hid_t) h5_root_id,dataset_id,group_id
integer ierr,i
integer(HSIZE_T), dimension(ndims) :: h_dims
character*100 errmsg


do i=1,ndims
  h_dims(i)=dims(i)
enddo

call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F,h5_root_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(read_real8_array) : h5fopen_f returned",I6)')ierr
  goto 100
endif
call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
if (ierr.ne.0) then
  write(errmsg,'("Error(read_real8_array) : h5gopen_f returned",I6)')ierr
  goto 100
endif
call h5dopen_f(group_id,trim(nm),dataset_id,ierr)
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

!---------------------------!
!    write_integer_array    !
!---------------------------!
subroutine write_integer_array(a,ndims,dims,fname,path,nm)
use hdf5
implicit none
integer, intent(in) :: ndims
integer, intent(in) :: dims(ndims)
integer, intent(in) :: a(*)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm


integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr,i
integer(HSIZE_T), dimension(ndims) :: h_dims

do i=1,ndims
  h_dims(i)=dims(i)
enddo
call h5fopen_f(trim(fname),H5F_ACC_RDWR_F,h5_root_id,ierr)
call h5screate_simple_f(ndims,h_dims,dataspace_id,ierr)
call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
call h5dcreate_f(group_id,trim(nm),H5T_NATIVE_INTEGER,dataspace_id, &
  dataset_id,ierr)
call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,a,h_dims,ierr)
call h5gclose_f(group_id,ierr)
call h5sclose_f(dataspace_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end

!--------------------------!
!    read_integer_array    !
!--------------------------!
subroutine read_integer_array(a,ndims,dims,fname,path,nm)
use hdf5
implicit none
integer, intent(in) :: ndims
integer, intent(in) :: dims(ndims)
integer, intent(out) :: a(*)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm


integer(hid_t) h5_root_id,dataset_id,group_id
integer ierr,i
integer(HSIZE_T), dimension(ndims) :: h_dims

do i=1,ndims
  h_dims(i)=dims(i)
enddo
call h5fopen_f(trim(fname),H5F_ACC_RDONLY_F,h5_root_id,ierr)
call h5gopen_f(h5_root_id,trim(path),group_id,ierr)
call h5dopen_f(group_id,trim(nm),dataset_id,ierr)
call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,a,h_dims,ierr)
call h5gclose_f(group_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end
#endif

#ifndef _MPI_
subroutine pstop
implicit none
stop
return
end
#endif


