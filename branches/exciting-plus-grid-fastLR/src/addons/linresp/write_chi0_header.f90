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
ie1=0
call h5fcreate_f(trim(fchi0),H5F_ACC_TRUNC_F,h5_root_id,ierr)
!call h5gcreate_f(h5_root_id,'parameters',h5_tmp_id,ierr)
!call h5gclose_f(h5_tmp_id,ierr)
!if (wannier_chi0_chi) then
!  call h5gcreate_f(h5_root_id,'wannier',h5_tmp_id,ierr)
!  call h5gclose_f(h5_tmp_id,ierr)
!endif
call h5gcreate_f(h5_root_id,'iw',h5_w_id,ierr)
do i=1,nepts
  write(c8,'(I8.8)')i
  call h5gcreate_f(h5_w_id,c8,h5_iw_id,ierr)
  call h5gclose_f(h5_iw_id,ierr)
enddo
call h5gclose_f(h5_w_id,ierr)
call h5fclose_f(h5_root_id,ierr)
!call write_integer(nepts,1,trim(fchi0),'/parameters','nepts')
!call write_integer(lr_igq0,1,trim(fchi0),'/parameters','lr_igq0')
!call write_integer(gshme1,1,trim(fchi0),'/parameters','gshme1')
!call write_integer(gshme2,1,trim(fchi0),'/parameters','gshme2')
!call write_integer(gvecme1,1,trim(fchi0),'/parameters','gvecme1')
!call write_integer(gvecme2,1,trim(fchi0),'/parameters','gvecme2')
!call write_integer(ngvecme,1,trim(fchi0),'/parameters','ngvecme')
!call write_real8(vq0l,3,trim(fchi0),'/parameters','vq0l')
!call write_real8(vq0rl,3,trim(fchi0),'/parameters','vq0rl')
!call write_real8(vq0c,3,trim(fchi0),'/parameters','vq0c')
!call write_real8(vq0rc,3,trim(fchi0),'/parameters','vq0rc')
!call write_integer(ie1,1,trim(fchi0),'/parameters','ie1')
return
end