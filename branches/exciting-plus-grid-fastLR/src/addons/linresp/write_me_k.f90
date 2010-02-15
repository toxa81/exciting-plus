subroutine write_me_k(ikloc,fname,me)
#ifdef _HDF5_
use modmain
implicit none
integer, intent(in) :: ikloc
character*(*), intent(in) :: fname 
complex(8), intent(in) :: me(ngvecme,nmegqblhmax)

character*100 path
integer ik
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
write(path,'("/kpoints/",I8.8)')ik
call write_integer(idxkq(1,ik),1,trim(fname),trim(path),'kq')
call write_integer(nmegqblh(ikloc),1,trim(fname),trim(path),'nmegqblh')
if (nmegqblh(ikloc).gt.0) then
  call write_integer_array(bmegqblh(1,1,ikloc),2,(/2,nmegqblh(ikloc)/), &
    trim(fname),trim(path),'bmegqblh')
  call write_real8_array(me,3,(/2,ngvecme,nmegqblh(ikloc)/), &
    trim(fname),trim(path),'megqblh')
endif
if (wannier_megq) then
  call write_real8_array(wann_c(1,1,ikloc),3,(/2,nwann,nstsv/), &
    trim(fname),trim(path),'wann_c_k')
  call write_real8_array(wann_c(1,1,ikloc+nkptnrloc),3,(/2,nwann,nstsv/), &
    trim(fname),trim(path),'wann_c_kq')
endif   
#endif
return
end
