#ifdef _HDF5_
subroutine readmegqblh(qnm)
use modmain
use mod_nrkp
implicit none
character(*), intent(in) :: qnm
integer i,j,ik,ikloc
integer nkptnr_,nspinor_
character*100 fme,fmek,path

if (allocated(evalsvnr)) deallocate(evalsvnr)
allocate(evalsvnr(nstsv,nkptnr))
if (allocated(occsvnr)) deallocate(occsvnr)
allocate(occsvnr(nstsv,nkptnr))

fme=trim(qnm)//"_me.hdf5"
if (mpi_grid_root((/dim_k,dim_b/))) then
  call read_integer(nkptnr_,1,trim(fme),'/parameters','nkptnr')
  call read_integer(nmegqblhmax,1,trim(fme),'/parameters','nmegqblhmax')
  call read_integer(lr_igq0,1,trim(fme),'/parameters','lr_igq0')
  call read_integer(gshme1,1,trim(fme),'/parameters','gshme1')
  call read_integer(gshme2,1,trim(fme),'/parameters','gshme2')
  call read_integer(gvecme1,1,trim(fme),'/parameters','gvecme1')
  call read_integer(gvecme2,1,trim(fme),'/parameters','gvecme2')
  call read_integer(ngvecme,1,trim(fme),'/parameters','ngvecme')
  call read_integer(nspinor_,1,trim(fme),'/parameters','nspinor')
  call read_real8(vq0l,3,trim(fme),'/parameters','vq0l')
  call read_real8(vq0rl,3,trim(fme),'/parameters','vq0rl')
  call read_real8(vq0c,3,trim(fme),'/parameters','vq0c')
  call read_real8(vq0rc,3,trim(fme),'/parameters','vq0rc')
  call read_real8_array(evalsvnr,2,(/nstsv,nkptnr/), &
      trim(fme),'/parameters','evalsvnr')
  call read_real8_array(occsvnr,2,(/nstsv,nkptnr/), &
      trim(fme),'/parameters','occsvnr')  
  if (nkptnr_.ne.nkptnr) then
    write(*,*)
    write(*,'("Error(readmegq): k-mesh was changed")')
    write(*,*)
    call pstop
  endif
  if (nspinor_.ne.nspinor) then
    write(*,*)
    write(*,'("Error(readmegq): number of spin components was changed")')
    write(*,*)
    call pstop
  endif    
endif
call mpi_grid_bcast(gshme1,dims=(/dim_k,dim_b/))
call mpi_grid_bcast(gshme2,dims=(/dim_k,dim_b/))
call mpi_grid_bcast(gvecme1,dims=(/dim_k,dim_b/))
call mpi_grid_bcast(gvecme2,dims=(/dim_k,dim_b/))
call mpi_grid_bcast(ngvecme,dims=(/dim_k,dim_b/))
call mpi_grid_bcast(nmegqblhmax,dims=(/dim_k,dim_b/))
call mpi_grid_bcast(lr_igq0,dims=(/dim_k,dim_b/))
call mpi_grid_bcast(vq0l(1),3,dims=(/dim_k,dim_b/))
call mpi_grid_bcast(vq0rl(1),3,dims=(/dim_k,dim_b/))
call mpi_grid_bcast(vq0c(1),3,dims=(/dim_k,dim_b/))
call mpi_grid_bcast(vq0rc(1),3,dims=(/dim_k,dim_b/))
call mpi_grid_bcast(evalsvnr(1,1),nstsv*nkptnr,dims=(/dim_k,dim_b/))
call mpi_grid_bcast(occsvnr(1,1),nstsv*nkptnr,dims=(/dim_k,dim_b/))

allocate(idxkq(1,nkptnr))
idxkq=0
allocate(nmegqblh(nkptnrloc))
allocate(bmegqblh(2,nmegqblhmax,nkptnrloc))
allocate(megqblh(ngvecme,nmegqblhmax,nkptnrloc))

! read matrix elements
if (.not.split_megq_file) then
  if (mpi_grid_side(dims=(/dim_k,dim_q/))) then
    do i=0,mpi_grid_size(1)-1
    do j=0,mpi_grid_size(3)-1
      if (mpi_grid_x(1).eq.i.and.mpi_grid_x(3).eq.j) then
        do ikloc=1,nkptnrloc
          ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
          write(path,'("/kpoints/",I8.8)')ik         
          call read_integer(idxkq(1,ik),1,trim(fme),trim(path),'kq')
          call read_integer(nmegqblh(ikloc),1,trim(fme),trim(path),'nmegqblh')
          if (nmegqblh(ikloc).gt.0) then
            call read_integer_array(bmegqblh(1,1,ikloc),2,(/2,nmegqblh(ikloc)/), &
              trim(fme),trim(path),'bmegqblh')
            call read_real8_array(megqblh(1,1,ikloc),3,(/2,ngvecme,nmegqblh(ikloc)/), &
              trim(fme),trim(path),'megqblh')
          endif
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
      call read_integer(idxkq(1,ik),1,trim(fmek),trim(path),'kq')
      call read_integer(nmegqblh(ikloc),1,trim(fmek),trim(path),'nmegqblh')
      if (nmegqblh(ikloc).gt.0) then
        call read_integer_array(bmegqblh(1,1,ikloc),2,(/2,nmegqblh(ikloc)/), &
          trim(fmek),trim(path),'bmegqblh')
        call read_real8_array(megqblh(1,1,ikloc),3,(/2,ngvecme,nmegqblh(ikloc)/), &
          trim(fmek),trim(path),'megqblh')
      endif
    enddo
  endif
endif
call mpi_grid_barrier(dims=(/dim_k,dim_b/))
call mpi_grid_reduce(idxkq(1,1),nkptnr,dims=(/dim_k/),all=.true.)
call mpi_grid_bcast(idxkq(1,1),nkptnr,dims=(/dim_b/))
call mpi_grid_bcast(nmegqblh(1),nkptnrloc,dims=(/dim_b/))
call mpi_grid_bcast(bmegqblh(1,1,1),2*nmegqblhmax*nkptnrloc,dims=(/dim_b/))
call mpi_grid_bcast(megqblh(1,1,1),ngvecme*nmegqblhmax*nkptnrloc,dims=(/dim_b/))
return
end
#endif