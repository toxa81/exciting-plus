subroutine readmegqwan(qnm)
use modmain
implicit none
character(*), intent(in) :: qnm
integer i,j,ik,ikloc
integer nkptnr_,nspinor_
character*100 fme,fmek,path

if (allocated(wann_occ)) deallocate(wann_occ)
allocate(wann_occ(nwann))
if (allocated(wann_c)) deallocate(wann_c)
allocate(wann_c(nwann,nstsv,2*nkptnrloc))

fme=trim(qnm)//"_me.hdf5"
if (mpi_grid_root((/dim_k,dim_b/))) then
  call read_real8(wann_occ,nwann,trim(fme),'/wannier','wann_occ')
  !call read_integer(ntrmegqwan,1,trim(fme),'/wannier','ntrmegqwan')
  call read_integer(nmegqwan,1,trim(fme),'/wannier','nmegqwan')
endif
call mpi_grid_bcast(wann_occ(1),nwann,dims=(/dim_k,dim_b/))
!call mpi_grid_bcast(ntrmegqwan,dims=(/dim_k,dim_b/))
call mpi_grid_bcast(nmegqwan,dims=(/dim_k,dim_b/))

!allocate(itrmegqwan(3,ntrmegqwan))
allocate(megqwan(nmegqwan,ngvecme))
!allocate(bmegqwan(2,nwann*nwann))
if (mpi_grid_root((/dim_k,dim2/))) then
  !call read_integer_array(itrmegqwan,2,(/3,ntrmegqwan/),trim(fme),'/wannier','itrmegqwan')
  call read_real8_array(megqwan,3,(/2,nmegqwan,ngvecme/), &
    trim(fme),'/wannier','megqwan')
  !call read_integer_array(bmegqwan,2,(/2,nmegqwan/),trim(fme),'/wannier','bmegqwan')
endif
!call mpi_grid_bcast(itrmegqwan(1,1),3*ntrmegqwan,dims=(/dim_k,dim2/))
!call mpi_grid_bcast(bmegqwan(1,1),2*nmegqwan,dims=(/dim_k,dim2/))
call mpi_grid_bcast(megqwan(1,1),nmegqwan*ngvecme,dims=(/dim_k,dim2/))

! read matrix elements
if (.not.split_megq_file) then
  if (mpi_grid_side(dims=(/dim_k,dim_q/))) then
    do i=0,mpi_grid_size(1)-1
    do j=0,mpi_grid_size(3)-1
      if (mpi_grid_x(1).eq.i.and.mpi_grid_x(3).eq.j) then
        do ikloc=1,nkptnrloc
          ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
          write(path,'("/kpoints/",I8.8)')ik         
          call read_real8_array(wann_c(1,1,ikloc),3,(/2,nwann,nstsv/), &
            trim(fme),trim(path),'wann_c_k')
          call read_real8_array(wann_c(1,1,ikloc+nkptnrloc),3,(/2,nwann,nstsv/), &
            trim(fme),trim(path),'wann_c_kq')
!          if (lwannopt) then
!            call read_real8_array(pmat(1,1,1,ikloc),4,(/2,3,nstsv,nstsv/), &
!              trim(fme),trim(path),'pmat')          
!          endif
        enddo
      endif
      if (.not.parallel_read) call mpi_grid_barrier(dims=(/dim_k,dim_q/))
    enddo !j
    enddo !i
  endif
else
  if (mpi_grid_root((/dim_b/))) then
    do ikloc=1,nkptnrloc
      ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
      write(fmek,'("_me_k_",I8.8)')ik
      fmek=trim(qnm)//trim(fmek)//".hdf5"
      write(path,'("/kpoints/",I8.8)')ik
      call read_real8_array(wann_c(1,1,ikloc),3,(/2,nwann,nstsv/), &
        trim(fmek),trim(path),'wann_c_k')
      call read_real8_array(wann_c(1,1,ikloc+nkptnrloc),3,(/2,nwann,nstsv/), &
        trim(fmek),trim(path),'wann_c_kq')
!      if (lwannopt) then
!        call read_real8_array(pmat(1,1,1,ikloc),4,(/2,3,nstsv,nstsv/), &
!          trim(fmek),trim(path),'pmat')          
!      endif
    enddo
  endif
endif
call mpi_grid_barrier(dims=(/dim_k,dim_b/))
call mpi_grid_bcast(wann_c(1,1,1),nwann*nstsv*2*nkptnrloc,dims=(/dim2/))
!if (lwannopt) then
!  call d_bcast_cart(comm_cart_010,pmat,2*3*nstsv*nstsv*nkptnr_loc)
!endif
return
end