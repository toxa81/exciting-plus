subroutine write_me_header(qnm)
use modmain
use hdf5
implicit none
character(*), intent(in) :: qnm
integer ik,ikloc,ierr,complete
integer(hid_t) h5_root_id
integer(hid_t) h5_kpoints_id
integer(hid_t) h5_kpoint_id
integer(hid_t) h5_tmp_id
character*8 c8
character*100 fme,fmek

fme=trim(qnm)//"_me.hdf5"

if (mpi_grid_root(dims=(/dim_k,dim2/))) then
  call h5fcreate_f(trim(fme),H5F_ACC_TRUNC_F,h5_root_id,ierr)
  call h5gcreate_f(h5_root_id,'parameters',h5_tmp_id,ierr)
  call h5gclose_f(h5_tmp_id,ierr)
  if (wannier_megq) then
    call h5gcreate_f(h5_root_id,'wannier',h5_tmp_id,ierr)
    call h5gclose_f(h5_tmp_id,ierr)
  endif
  call h5gcreate_f(h5_root_id,'kpoints',h5_kpoints_id,ierr)
  do ik=1,nkptnr
    write(c8,'(I8.8)')ik
    call h5gcreate_f(h5_kpoints_id,c8,h5_kpoint_id,ierr)
    call h5gclose_f(h5_kpoint_id,ierr)
  enddo
  call h5gclose_f(h5_kpoints_id,ierr)
  call h5fclose_f(h5_root_id,ierr)
  call write_integer(nkptnr,1,trim(fme),'/parameters','nkptnr')
  call write_integer(nmegqblhmax,1,trim(fme),'/parameters','nmegqblhmax')
  call write_integer(lr_igq0,1,trim(fme),'/parameters','lr_igq0')
  call write_integer(gshme1,1,trim(fme),'/parameters','gshme1')
  call write_integer(gshme2,1,trim(fme),'/parameters','gshme2')
  call write_integer(gvecme1,1,trim(fme),'/parameters','gvecme1')
  call write_integer(gvecme2,1,trim(fme),'/parameters','gvecme2')
  call write_integer(ngvecme,1,trim(fme),'/parameters','ngvecme')
  call write_integer(nspinor,1,trim(fme),'/parameters','nspinor')
  call write_real8(vq0l,3,trim(fme),'/parameters','vq0l')
  call write_real8(vq0rl,3,trim(fme),'/parameters','vq0rl')
  call write_real8(vq0c,3,trim(fme),'/parameters','vq0c')
  call write_real8(vq0rc,3,trim(fme),'/parameters','vq0rc')
  call write_real8_array(lr_evalsvnr,2,(/nstsv,nkptnr/), &
      trim(fme),'/parameters','evalsvnr')
  call write_real8_array(lr_occsvnr,2,(/nstsv,nkptnr/), &
      trim(fme),'/parameters','occsvnr')  
  complete=0
  call write_integer(complete,1,trim(fme),'/parameters','complete')
  if (wannier_megq) then
    call write_real8(wann_occ,nwann,trim(fme),'/wannier','wann_occ')
  endif
endif
if (mpi_grid_root(dims=(/dim2/)).and.split_megq_file) then
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    write(c8,'(I8.8)')ik
    fmek=trim(qnm)//"_me_k_"//c8//".hdf5"
    call h5fcreate_f(trim(fmek),H5F_ACC_TRUNC_F,h5_root_id,ierr)
    call h5gcreate_f(h5_root_id,'kpoints',h5_kpoints_id,ierr)
    call h5gcreate_f(h5_kpoints_id,c8,h5_kpoint_id,ierr)
    call h5gclose_f(h5_kpoint_id,ierr)
    call h5gclose_f(h5_kpoints_id,ierr)
    call h5fclose_f(h5_root_id,ierr)
  enddo
endif

return
end