#ifdef _HDF5_
subroutine write_real8(a,n,fname,path,nm)
use hdf5
implicit none
integer, intent(in) :: n
real(8), intent(in) :: a(n)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm


integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr
integer(HSIZE_T), dimension(1) :: dims

dims(1)=n
call h5fopen_f(fname,H5F_ACC_RDWR_F,h5_root_id,ierr)
call h5screate_simple_f(1,dims,dataspace_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dcreate_f(group_id,nm,H5T_NATIVE_DOUBLE,dataspace_id, &
  dataset_id,ierr)
call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,a,dims,ierr)
call h5gclose_f(group_id,ierr)
call h5sclose_f(dataspace_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end

subroutine read_real8(a,n,fname,path,nm)
use hdf5
implicit none
integer, intent(in) :: n
real(8), intent(out) :: a(n)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm

integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr
integer(HSIZE_T), dimension(1) :: dims

dims(1)=n
call h5fopen_f(fname,H5F_ACC_RDONLY_F,h5_root_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dopen_f(group_id,nm,dataset_id,ierr)
call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,a,dims,ierr)
call h5gclose_f(group_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end


subroutine write_real8_array(a,ndims,dims,fname,path,nm)
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
integer(HSIZE_T), dimension(ndims) :: h_dims

do i=1,ndims
  h_dims(i)=dims(i)
enddo
call h5fopen_f(fname,H5F_ACC_RDWR_F,h5_root_id,ierr)
call h5screate_simple_f(ndims,h_dims,dataspace_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dcreate_f(group_id,nm,H5T_NATIVE_DOUBLE,dataspace_id, &
  dataset_id,ierr)
call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,a,h_dims,ierr)
call h5gclose_f(group_id,ierr)
call h5sclose_f(dataspace_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end

subroutine read_real8_array(a,ndims,dims,fname,path,nm)
use hdf5
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

do i=1,ndims
  h_dims(i)=dims(i)
enddo

call h5fopen_f(fname,H5F_ACC_RDONLY_F,h5_root_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dopen_f(group_id,nm,dataset_id,ierr)
call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,a,h_dims,ierr)
call h5gclose_f(group_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end


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
call h5fopen_f(fname,H5F_ACC_RDWR_F,h5_root_id,ierr)
call h5screate_simple_f(ndims,h_dims,dataspace_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dcreate_f(group_id,nm,H5T_NATIVE_INTEGER,dataspace_id, &
  dataset_id,ierr)
call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,a,h_dims,ierr)
call h5gclose_f(group_id,ierr)
call h5sclose_f(dataspace_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end


subroutine read_integer_array(a,ndims,dims,fname,path,nm)
use hdf5
implicit none
integer, intent(in) :: ndims
integer, intent(in) :: dims(ndims)
integer, intent(out) :: a(*)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm


integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr,i
integer(HSIZE_T), dimension(ndims) :: h_dims

do i=1,ndims
  h_dims(i)=dims(i)
enddo
call h5fopen_f(fname,H5F_ACC_RDONLY_F,h5_root_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dopen_f(group_id,nm,dataset_id,ierr)
call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,a,h_dims,ierr)
call h5gclose_f(group_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end

subroutine write_integer(a,n,fname,path,nm)
use hdf5
implicit none
integer, intent(in) :: n
integer, intent(in) :: a(n)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm


integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr
integer(HSIZE_T), dimension(1) :: dims

dims(1)=n
call h5fopen_f(fname,H5F_ACC_RDWR_F,h5_root_id,ierr)
call h5screate_simple_f(1,dims,dataspace_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dcreate_f(group_id,nm,H5T_NATIVE_INTEGER,dataspace_id, &
  dataset_id,ierr)
call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,a,dims,ierr)
call h5gclose_f(group_id,ierr)
call h5sclose_f(dataspace_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end

subroutine rewrite_integer(a,n,fname,path,nm)
use hdf5
implicit none
integer, intent(in) :: n
integer, intent(in) :: a(n)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm


integer(hid_t) h5_root_id,dataset_id,group_id
integer ierr
integer(HSIZE_T), dimension(1) :: dims

dims(1)=n
call h5fopen_f(fname,H5F_ACC_RDWR_F,h5_root_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dopen_f(group_id,nm,dataset_id,ierr)
call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,a,dims,ierr)
call h5gclose_f(group_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end

subroutine read_integer(a,n,fname,path,nm)
use hdf5
implicit none
integer, intent(in) :: n
integer, intent(out) :: a(n)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm

integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr
integer(HSIZE_T), dimension(1) :: dims

dims(1)=n
call h5fopen_f(fname,H5F_ACC_RDONLY_F,h5_root_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dopen_f(group_id,nm,dataset_id,ierr)
call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,a,dims,ierr)
call h5gclose_f(group_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end
#endif