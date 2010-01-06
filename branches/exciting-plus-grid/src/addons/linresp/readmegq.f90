subroutine readmegq(qnm)
use modmain
implicit none
character(*), intent(in) :: qnm
integer i,j,ik,ikloc
integer nkptnr_,nspinor_
character*100 fme,fmek,path

if (allocated(lr_evalsvnr)) deallocate(lr_evalsvnr)
allocate(lr_evalsvnr(nstsv,nkptnr))
if (allocated(lr_occsvnr)) deallocate(lr_occsvnr)
allocate(lr_occsvnr(nstsv,nkptnr))
if (allocated(wann_occ)) deallocate(wann_occ)
allocate(wann_occ(nwann))

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
  if (wannier_megq) then
    call read_real8(wann_occ,nwann,trim(fme),'/wannier','wann_occ')
    call read_integer(ntrmegqwan,1,trim(fme),'/wannier','ntrmegqwan')
    call read_integer(nmegqwan,1,trim(fme),'/wannier','nmegqwan')
  endif 
  call read_real8_array(lr_evalsvnr,2,(/nstsv,nkptnr/), &
      trim(fme),'/parameters','evalsvnr')
  call read_real8_array(lr_occsvnr,2,(/nstsv,nkptnr/), &
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
if (wannier_megq) then
  call mpi_grid_bcast(wann_occ(1),nwann,dims=(/dim_k,dim_b/))
  call mpi_grid_bcast(ntrmegqwan,dims=(/dim_k,dim_b/))
  call mpi_grid_bcast(nmegqwan,dims=(/dim_k,dim_b/))
endif  
call mpi_grid_bcast(lr_evalsvnr(1,1),nstsv*nkptnr,dims=(/dim_k,dim_b/))
call mpi_grid_bcast(lr_occsvnr(1,1),nstsv*nkptnr,dims=(/dim_k,dim_b/))

if (wannier_megq) then
  allocate(itrmegqwan(3,ntrmegqwan))
  allocate(megqwan(nmegqwan,ntrmegqwan,ngvecme))
  allocate(bmegqwan(2,nwann*nwann))
  if (mpi_grid_root((/dim_k,dim2/))) then
    call read_integer_array(itrmegqwan,2,(/3,ntrmegqwan/),trim(fme),'/wannier','itrmegqwan')
    call read_real8_array(megqwan,4,(/2,nmegqwan,ntrmegqwan,ngvecme/), &
      trim(fme),'/wannier','megqwan')
    call read_integer_array(bmegqwan,2,(/2,nmegqwan/),trim(fme),'/wannier','bmegqwan')
  endif
  call mpi_grid_bcast(itrmegqwan(1,1),3*ntrmegqwan,dims=(/dim_k,dim2/))
  call mpi_grid_bcast(bmegqwan(1,1),2*nmegqwan,dims=(/dim_k,dim2/))
  call mpi_grid_bcast(megqwan(1,1,1),nmegqwan*ntrmegqwan*ngvecme,dims=(/dim_k,dim2/))
endif

allocate(idxkq(1,nkptnr))
idxkq=0
allocate(nmegqblh(nkptnrloc))
allocate(bmegqblh(2,nmegqblhmax,nkptnrloc))
allocate(megqblh(ngvecme,nmegqblhmax,nkptnrloc))
if (wannier_megq) then
  if (allocated(wann_c)) deallocate(wann_c)
  allocate(wann_c(nwann,nstsv,2*nkptnrloc))
endif  

if (wproc) then
  write(150,*)
  if (nspinor_.eq.2) then
    write(150,'("matrix elements were calculated for spin-polarized case")')
  endif
  write(150,'("matrix elements were calculated for: ")')
  write(150,'("  G-shells  : ",I4," to ", I4)')gshme1,gshme2
  write(150,'("  G-vectors : ",I4," to ", I4)')gvecme1,gvecme2
  write(150,'("Reading matrix elements")')
  call flushifc(150)
endif

! read matrix elements
call timer_start(1,reset=.true.)
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
          if (wannier_megq) then
            call read_real8_array(wann_c(1,1,ikloc),3,(/2,nwann,nstsv/), &
              trim(fme),trim(path),'wann_c_k')
            call read_real8_array(wann_c(1,1,ikloc+nkptnrloc),3,(/2,nwann,nstsv/), &
              trim(fme),trim(path),'wann_c_kq')
          endif 
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
      call read_integer(idxkq(1,ik),1,trim(fmek),trim(path),'kq')
      call read_integer(nmegqblh(ikloc),1,trim(fmek),trim(path),'nmegqblh')
      if (nmegqblh(ikloc).gt.0) then
        call read_integer_array(bmegqblh(1,1,ikloc),2,(/2,nmegqblh(ikloc)/), &
          trim(fmek),trim(path),'bmegqblh')
        call read_real8_array(megqblh(1,1,ikloc),3,(/2,ngvecme,nmegqblh(ikloc)/), &
          trim(fmek),trim(path),'megqblh')
      endif
      if (wannier_megq) then
        call read_real8_array(wann_c(1,1,ikloc),3,(/2,nwann,nstsv/), &
          trim(fmek),trim(path),'wann_c_k')
        call read_real8_array(wann_c(1,1,ikloc+nkptnrloc),3,(/2,nwann,nstsv/), &
          trim(fmek),trim(path),'wann_c_kq')
      endif    
!      if (lwannopt) then
!        call read_real8_array(pmat(1,1,1,ikloc),4,(/2,3,nstsv,nstsv/), &
!          trim(fmek),trim(path),'pmat')          
!      endif
    enddo
  endif
endif
call mpi_grid_barrier(dims=(/dim_k,dim_b/))
call timer_stop(1)
if (wproc) then
   write(150,'("Done in ",F8.2," seconds")')timer_get_value(1)
  call flushifc(150)
endif
call mpi_grid_reduce(idxkq(1,1),nkptnr,dims=(/dim_k/),all=.true.)
call mpi_grid_bcast(idxkq(1,1),nkptnr,dims=(/dim_b/))
call mpi_grid_bcast(nmegqblh(1),nkptnrloc,dims=(/dim_b/))
call mpi_grid_bcast(bmegqblh(1,1,1),2*nmegqblhmax*nkptnrloc,dims=(/dim_b/))
call mpi_grid_bcast(megqblh(1,1,1),ngvecme*nmegqblhmax*nkptnrloc,dims=(/dim_b/))
if (wannier_megq) then
  call mpi_grid_bcast(wann_c(1,1,1),nwann*nstsv*2*nkptnrloc,dims=(/dim2/))
endif
!if (lwannopt) then
!  call d_bcast_cart(comm_cart_010,pmat,2*3*nstsv*nstsv*nkptnr_loc)
!endif
return
end