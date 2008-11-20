subroutine init3
#ifdef _HDF5_
use hdf5
use modmain
implicit none
integer ierr
integer(hid_t) h5_root_id
integer(hid_t) h5_kpoints_id
integer(hid_t) h5_kpoint_id
integer(hid_t) h5_tmp_id
character*4 c4

integer ik
call h5fcreate_f('exciting.hdf5',H5F_ACC_TRUNC_F,h5_root_id,ierr)
call h5gcreate_f(h5_root_id,'dimensions',h5_tmp_id,ierr)
call h5gclose_f(h5_tmp_id,ierr)
call h5gcreate_f(h5_root_id,'kpoints',h5_kpoints_id,ierr)
do ik=1,nkpt
  write(c4,'(I4.4)')ik
  call h5gcreate_f(h5_kpoints_id,c4,h5_kpoint_id,ierr)
  call h5gclose_f(h5_kpoint_id,ierr)
enddo
call h5gclose_f(h5_kpoints_id,ierr)
call h5fclose_f(h5_root_id,ierr)
#endif  
return
end

subroutine write_integer(a,n,path,nm)
#ifdef _HDF5_
use hdf5
implicit none
integer, intent(in) :: n
integer, intent(in) :: a(n)
character(*), intent(in) :: path
character(*), intent(in) :: nm


integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr
integer(HSIZE_T), dimension(1) :: dims

dims(1)=n
call h5fopen_f('exciting.hdf5',H5F_ACC_RDWR_F,h5_root_id,ierr)
call h5screate_simple_f(1,dims,dataspace_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dcreate_f(group_id,nm,H5T_NATIVE_INTEGER,dataspace_id, &
  dataset_id,ierr)
call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,a,dims,ierr)
call h5gclose_f(group_id,ierr)
call h5sclose_f(dataspace_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
#else
write(*,'("Error(hdf5_io): not compiled with HDF5 support")')
call pstop
#endif
return
end

subroutine read_integer(a,n,path,nm)
#ifdef _HDF5_
use hdf5
implicit none
integer, intent(in) :: n
integer, intent(out) :: a(n)
character(*), intent(in) :: path
character(*), intent(in) :: nm


integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr
integer(HSIZE_T), dimension(1) :: dims

dims(1)=n
call h5fopen_f('exciting.hdf5',H5F_ACC_RDWR_F,h5_root_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dopen_f(group_id,nm,dataset_id,ierr)
call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,a,dims,ierr)
call h5gclose_f(group_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
#else
write(*,'("Error(hdf5_io): not compiled with HDF5 support")')
call pstop
#endif
return
end

subroutine write_real8(a,n,path,nm)
#ifdef _HDF5_
use hdf5
implicit none
integer, intent(in) :: n
real(8), intent(in) :: a(n)
character(*), intent(in) :: path
character(*), intent(in) :: nm


integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr
integer(HSIZE_T), dimension(1) :: dims

dims(1)=n
call h5fopen_f('exciting.hdf5',H5F_ACC_RDWR_F,h5_root_id,ierr)
call h5screate_simple_f(1,dims,dataspace_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dcreate_f(group_id,nm,H5T_NATIVE_DOUBLE,dataspace_id, &
  dataset_id,ierr)
call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,a,dims,ierr)
call h5gclose_f(group_id,ierr)
call h5sclose_f(dataspace_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
#else
write(*,'("Error(hdf5_io): not compiled with HDF5 support")')
call pstop
#endif
return
end

subroutine write_complex16(a,n,path,nm)
#ifdef _HDF5_
implicit none
integer, intent(in) :: n
complex(8), intent(in) :: a(n)
character(*), intent(in) :: path
character(*), intent(in) :: nm

call write_real8(a,n*2,path,nm)
#else
write(*,'("Error(hdf5_io): not compiled with HDF5 support")')
call pstop
#endif

return
end

subroutine read_real8(a,n,path,nm)
#ifdef _HDF5_
use hdf5
implicit none
integer, intent(in) :: n
real(8), intent(out) :: a(n)
character(*), intent(in) :: path
character(*), intent(in) :: nm


integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr
integer(HSIZE_T), dimension(1) :: dims

dims(1)=n
call h5fopen_f('exciting.hdf5',H5F_ACC_RDWR_F,h5_root_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dopen_f(group_id,nm,dataset_id,ierr)
call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,a,dims,ierr)
call h5gclose_f(group_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
#else
write(*,'("Error(hdf5_io): not compiled with HDF5 support")')
call pstop
#endif

return
end

subroutine read_complex16(a,n,path,nm)
#ifdef _HDF5_
implicit none
integer, intent(in) :: n
complex(8), intent(out) :: a(n)
character(*), intent(in) :: path
character(*), intent(in) :: nm

call read_real8(a,n*2,path,nm)

#else
write(*,'("Error(hdf5_io): not compiled with HDF5 support")')
call pstop
#endif

return
end

