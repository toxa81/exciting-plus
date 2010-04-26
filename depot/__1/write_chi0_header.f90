#ifdef _HDF5_
subroutine write_chi0_header(qnm)
use modmain
implicit none
character(*), intent(in) :: qnm
character*100 fchi0
integer i
character*8 c8
fchi0=trim(qnm)//"_chi0.hdf5"
call hdf5_create_file(trim(fchi0))
call hdf5_create_group(trim(fchi0),'/','iw')
do i=1,lr_nw
  write(c8,'(I8.8)')i
  call hdf5_create_group(trim(fchi0),'/iw',c8)
enddo
return
end
#endif
