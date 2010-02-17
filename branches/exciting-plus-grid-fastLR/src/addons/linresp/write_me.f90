subroutine write_me(qnm)
use modmain
implicit none
character(*), intent(in) :: qnm
integer ikloc,ik,i
character*100 fme,fmek

fme=trim(qnm)//"_me.hdf5"

if (mpi_grid_root(dims=(/dim2/))) then
  if (.not.split_megq_file) then
    do i=0,mpi_grid_size(dim_k)-1
      if (mpi_grid_x(dim_k).eq.i) then
        do ikloc=1,nkptnrloc
          call write_me_k(ikloc,fme,megqblh(1,1,ikloc))
        enddo
      endif
      call mpi_grid_barrier(dims=(/dim_k/))
    enddo !i
  else
    do ikloc=1,nkptnrloc
      ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
      write(fmek,'("_me_k_",I8.8)')ik
      fmek=trim(qnm)//trim(fmek)//".hdf5"
      call write_me_k(ikloc,fmek,megqblh(1,1,ikloc))
    enddo
  endif !.not.spit_megq_file
endif !mpi_grid_root

if (wannier_megq) then
  if (mpi_grid_root((/dim_k,dim2/))) then
    !call write_integer(ntrmegqwan,1,trim(fme),'/wannier','ntrmegqwan')
    call write_integer(nmegqwan,1,trim(fme),'/wannier','nmegqwan')
    call write_real8_array(megqwan,3,(/2,nmegqwan,ngvecme/), &
      trim(fme),'/wannier','megqwan')
!    call write_integer_array(itrmegqwan,2,(/3,ntrmegqwan/),trim(fme),'/wannier','itrmegqwan')
!    call write_integer_array(bmegqwan,2,(/2,nmegqwan/),trim(fme),'/wannier','bmegqwan')
  endif
endif

return
end