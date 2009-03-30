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

!subroutine read_real8_array2(a,ndims,dims,fname,path,nm,h5_root_id)
!use hdf5
!implicit none
!integer, intent(in) :: ndims
!integer, intent(in) :: dims(ndims)
!real(8), intent(out) :: a(*)
!character(*), intent(in) :: fname
!character(*), intent(in) :: path
!character(*), intent(in) :: nm
!integer(hid_t),intent(in) :: h5_root_id
!
!integer(hid_t) dataspace_id,dataset_id,group_id
!integer ierr,i
!integer(HSIZE_T), dimension(ndims) :: h_dims
!
!do i=1,ndims
!  h_dims(i)=dims(i)
!enddo
!
!call h5gopen_f(h5_root_id,path,group_id,ierr)
!call h5dopen_f(group_id,nm,dataset_id,ierr)
!call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,a,h_dims,ierr)
!call h5gclose_f(group_id,ierr)
!call h5dclose_f(dataset_id,ierr)
!return
!end
!
!subroutine read_real8_array_p(a,ndims,dims,fname,path,nm,comm)
!use hdf5
!use mpi
!implicit none
!integer, intent(in) :: ndims
!integer, intent(in) :: dims(ndims)
!real(8), intent(out) :: a(*)
!character(*), intent(in) :: fname
!character(*), intent(in) :: path
!character(*), intent(in) :: nm
!integer, intent(in) :: comm
!
!integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id,plist_id
!integer ierr,i,info
!integer(HSIZE_T), dimension(ndims) :: h_dims
!
!info=MPI_INFO_NULL
!
!do i=1,ndims
!  h_dims(i)=dims(i)
!enddo
!
!call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,ierr)
!call h5pset_fapl_mpio_f(plist_id,comm,info,ierr)
!call h5fopen_f(fname,H5F_ACC_RDONLY_F,h5_root_id,ierr,access_prp=plist_id)
!call h5pclose_f(plist_id,ierr)
!
!call h5gopen_f(h5_root_id,path,group_id,ierr)
!call h5dopen_f(group_id,nm,dataset_id,ierr)
!
!call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,ierr)
!call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,ierr)
!!call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_INDEPENDENT_F,ierr)
!call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,a,h_dims,ierr,xfer_prp=plist_id)
!call h5pclose_f(plist_id,ierr)
!
!call h5gclose_f(group_id,ierr)
!call h5dclose_f(dataset_id,ierr)
!call h5fclose_f(h5_root_id,ierr)
!return
!end


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

!subroutine read_integer2(a,n,fname,path,nm,h5_root_id)
!use hdf5
!implicit none
!integer, intent(in) :: n
!integer, intent(out) :: a(n)
!character(*), intent(in) :: fname
!character(*), intent(in) :: path
!character(*), intent(in) :: nm
!integer(hid_t),intent(in):: h5_root_id
!
!integer(hid_t) dataspace_id,dataset_id,group_id
!integer ierr
!integer(HSIZE_T), dimension(1) :: dims
!
!dims(1)=n
!call h5gopen_f(h5_root_id,path,group_id,ierr)
!call h5dopen_f(group_id,nm,dataset_id,ierr)
!call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,a,dims,ierr)
!call h5gclose_f(group_id,ierr)
!call h5dclose_f(dataset_id,ierr)
!return
!end
#endif