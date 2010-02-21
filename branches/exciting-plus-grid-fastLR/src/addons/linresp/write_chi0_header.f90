subroutine write_chi0_header(qnm)
use modmain
use hdf5
implicit none
character(*), intent(in) :: qnm
integer ie1
integer(hid_t) h5_root_id
integer(hid_t) h5_w_id
integer(hid_t) h5_iw_id
integer(hid_t) h5_tmp_id
character*100 fchi0
integer i,ierr
character*8 c8

fchi0=trim(qnm)//"_chi0.hdf5"
call h5fcreate_f(trim(fchi0),H5F_ACC_TRUNC_F,h5_root_id,ierr)
if (ierr.ne.0) then 
  write(*,*)
  write(*,'("Error(write_chi0_header) : h5fcreate_f returned ",I4)')ierr
  write(*,'("  mpi_grid_x : ",5I4)')mpi_grid_x
  call pstop
endif
call h5gcreate_f(h5_root_id,'iw',h5_w_id,ierr)
if (ierr.ne.0) then 
  write(*,*)
  write(*,'("Error(write_chi0_header) : h5gcreate_f returned ",I4)')ierr
  write(*,'("  mpi_grid_x : ",5I4)')mpi_grid_x
  call pstop
endif
do i=1,nepts
  write(c8,'(I8.8)')i
  call h5gcreate_f(h5_w_id,c8,h5_iw_id,ierr)
  if (ierr.ne.0) then 
    write(*,*)
    write(*,'("Error(write_chi0_header) : h5gcreate_f returned ",I4)')ierr
    write(*,'("  mpi_grid_x : ",5I4)')mpi_grid_x
    write(*,'("  iw : ",I4)')i
    call pstop
  endif  
  call h5gclose_f(h5_iw_id,ierr)
  if (ierr.ne.0) then 
    write(*,*)
    write(*,'("Error(write_chi0_header) : h5gclose_f returned ",I4)')ierr
    write(*,'("  mpi_grid_x : ",5I4)')mpi_grid_x
    write(*,'("  iw : ",I4)')i
    call pstop
  endif  
enddo
call h5gclose_f(h5_w_id,ierr)
if (ierr.ne.0) then 
  write(*,*)
  write(*,'("Error(write_chi0_header) : h5gclose_f returned ",I4)')ierr
  write(*,'("  mpi_grid_x : ",5I4)')mpi_grid_x
  call pstop
endif
call h5fclose_f(h5_root_id,ierr)
if (ierr.ne.0) then 
  write(*,*)
  write(*,'("Error(write_chi0_header) : h5fclose_f returned ",I4)')ierr
  write(*,'("  mpi_grid_x : ",5I4)')mpi_grid_x
  call pstop
endif
return
end
