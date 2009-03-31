#ifdef _HDF5_
subroutine response_chi0(ivq0m,evalsvnr)
use modmain
use hdf5
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
integer, intent(in) :: ivq0m(3)
real(8), intent(in) :: evalsvnr(nstsv,nkptnr)
! local variables
complex(8), allocatable :: chi0_loc(:,:,:)
complex(8), allocatable :: me(:,:,:)
complex(4), allocatable :: c_chi0_loc(:,:,:)
complex(4), allocatable :: c_me(:,:,:)
integer i,j,ik,ie,nkptnr_,i1,i2,ikloc,ig1,ig2,nspinor_,ispn
integer idx0,bs
complex(8) wt
complex(4) c_wt
character*100 fname,path,qnm
character*8 c8
integer ierr
logical exist
integer ie1
logical, parameter :: lcmplx=.true.

! HDF5
integer(hid_t) h5_root_id
integer(hid_t) h5_w_id
integer(hid_t) h5_iw_id
integer(hid_t) h5_tmp_id

logical, external :: root_cart

if (wproc) then
  write(150,*)
  write(150,'("Calculation of KS polarisability chi0")')
  write(150,*)
  write(150,'("Energy mesh parameters:")')
  write(150,'("  maximum energy [eV] : ", F9.4)')maxomega
  write(150,'("  energy step    [eV] : ", F9.4)')domega
  write(150,'("  eta            [eV] : ", F9.4)')lr_eta
  call flushifc(150)
endif
  
! setup energy mesh
nepts=1+maxomega/domega
allocate(lr_w(nepts))
do i=1,nepts
  lr_w(i)=dcmplx(domega*(i-1),lr_eta)/ha2ev
enddo

call qname(ivq0m,qnm)
qnm="./"//trim(qnm)//"/"//trim(qnm)
fname=trim(qnm)//"_me.hdf5"

if (root_cart((/1,1,0/))) then
  call read_integer(nkptnr_,1,trim(fname),'/parameters','nkptnr')
  call read_integer(nmemax,1,trim(fname),'/parameters','nmemax')
  call read_integer(lr_igq0,1,trim(fname),'/parameters','lr_igq0')
  call read_integer(gshme1,1,trim(fname),'/parameters','gshme1')
  call read_integer(gshme2,1,trim(fname),'/parameters','gshme2')
  call read_integer(gvecme1,1,trim(fname),'/parameters','gvecme1')
  call read_integer(gvecme2,1,trim(fname),'/parameters','gvecme2')
  call read_integer(ngvecme,1,trim(fname),'/parameters','ngvecme')
  call read_integer(nspinor_,1,trim(fname),'/parameters','nspinor')
  call read_integer(spin_me,1,trim(fname),'/parameters','spin_me')
  call read_real8(vq0l,3,trim(fname),'/parameters','vq0l')
  call read_real8(vq0rl,3,trim(fname),'/parameters','vq0rl')
  call read_real8(vq0c,3,trim(fname),'/parameters','vq0c')
  call read_real8(vq0rc,3,trim(fname),'/parameters','vq0rc')
  if (nkptnr_.ne.nkptnr) then
    write(*,*)
    write(*,'("Error(response_chi0): k-mesh was changed")')
    write(*,*)
    call pstop
  endif
  if (nspinor_.ne.nspinor) then
    write(*,*)
    write(*,'("Error(response_chi0): number of spin components was changed")')
    write(*,*)
    call pstop
  endif    
endif

if (wproc) then
  if (nspinor_.eq.2) then
    write(150,'("  matrix elements were calculated for spin-polarized case")')
    if (spin_me.eq.1) write(150,'("  file contains spin-up matix elements")')
    if (spin_me.eq.2) write(150,'("  file contains spin-dn matix elements")')
    if (spin_me.eq.3) write(150,'("  file contains matix elements for both spins")')
  endif
  write(150,'("matrix elements were calculated for: ")')
  write(150,'("  G-shells  : ",I4," to ", I4)')gshme1,gshme2
  write(150,'("  G-vectors : ",I4," to ", I4)')gvecme1,gvecme2
endif

call i_bcast_cart(comm_cart_110,gshme1,1)
call i_bcast_cart(comm_cart_110,gshme2,1)
call i_bcast_cart(comm_cart_110,gvecme1,1)
call i_bcast_cart(comm_cart_110,gvecme2,1)
call i_bcast_cart(comm_cart_110,ngvecme,1)
call i_bcast_cart(comm_cart_110,spin_me,1)
call i_bcast_cart(comm_cart_110,nmemax,1)
call i_bcast_cart(comm_cart_110,lr_igq0,1)
call d_bcast_cart(comm_cart_110,vq0l,3)
call d_bcast_cart(comm_cart_110,vq0rl,3)
call d_bcast_cart(comm_cart_110,vq0c,3)
call d_bcast_cart(comm_cart_110,vq0rc,3)

if (spin_me.eq.3) then
  nspin_chi0=2
else
  nspin_chi0=1
endif

allocate(idxkq(1,nkptnr_loc))
allocate(nme(nkptnr_loc))
allocate(ime(3,nmemax,nkptnr_loc))
allocate(docc(nmemax,nkptnr_loc))
allocate(me(ngvecme,nmemax,nkptnr_loc))
if (lcmplx) allocate(c_me(ngvecme,nmemax,nkptnr_loc))

if (wproc) then
  write(150,'("Reading matrix elements")')
  call flushifc(150)
endif

call timer_reset(1)
call timer_start(1)
! read matrix elements
if (lsfio) then
if (root_cart((/0,1,0/))) then
#ifndef _PIO_
  do i=0,mpi_dims(1)-1
  do j=0,mpi_dims(3)-1
    if (mpi_x(1).eq.i.and.mpi_x(3).eq.j) then
#endif
      do ikloc=1,nkptnr_loc
        ik=ikptnrloc(mpi_x(1),1)+ikloc-1
        write(path,'("/kpoints/",I8.8)')ik
        call read_integer(idxkq(1,ikloc),1,trim(fname),trim(path),'kq')
        call read_integer(nme(ikloc),1,trim(fname),trim(path),'nme')
        call read_integer_array(ime(1,1,ikloc),2,(/3,nme(ikloc)/), &
          trim(fname),trim(path),'ime')
        call read_real8(docc(1,ikloc),nme(ikloc),trim(fname), &
          trim(path),'docc')
        call read_real8_array(me(1,1,ikloc),3,(/2,ngvecme,nme(ikloc)/), &
          trim(fname),trim(path),'me')
      enddo
#ifndef _PIO_      
    endif
    call barrier(comm_cart_101)
  enddo
  enddo
#endif
endif
else
  if (root_cart((/0,1,0/))) then
    do ikloc=1,nkptnr_loc
      ik=ikptnrloc(mpi_x(1),1)+ikloc-1
      write(fname,'("_me_k_",I8.8)')ik
      fname=trim(qnm)//trim(fname)//".hdf5"
      write(path,'("/kpoints/",I8.8)')ik
      call read_integer(idxkq(1,ikloc),1,trim(fname),trim(path),'kq')
      call read_integer(nme(ikloc),1,trim(fname),trim(path),'nme')
      call read_integer_array(ime(1,1,ikloc),2,(/3,nme(ikloc)/), &
        trim(fname),trim(path),'ime')
      call read_real8(docc(1,ikloc),nme(ikloc),trim(fname), &
        trim(path),'docc')
      call read_real8_array(me(1,1,ikloc),3,(/2,ngvecme,nme(ikloc)/), &
        trim(fname),trim(path),'me')
    enddo
  endif
endif

call barrier(comm_cart)
call timer_stop(1)
if (wproc) then
   write(150,'("Done in ",F8.2," seconds")')timer(1,2)
  call flushifc(150)
endif

call i_bcast_cart(comm_cart_010,idxkq,nkptnr_loc)
call i_bcast_cart(comm_cart_010,nme,nkptnr_loc)
call i_bcast_cart(comm_cart_010,ime,3*nmemax*nkptnr_loc)
call d_bcast_cart(comm_cart_010,docc,nmemax*nkptnr_loc)
call d_bcast_cart(comm_cart_010,me,2*ngvecme*nmemax*nkptnr_loc)

if (lcmplx) c_me=cmplx(me)

ie1=0
fname=trim(qnm)//"_chi0.hdf5"
if (root_cart((/1,1,0/))) then
  inquire(file=trim(fname),exist=exist)
  if (.not.exist) then
    call h5fcreate_f(trim(fname),H5F_ACC_TRUNC_F,h5_root_id,ierr)
    call h5gcreate_f(h5_root_id,'parameters',h5_tmp_id,ierr)
    call h5gclose_f(h5_tmp_id,ierr)
    call h5gcreate_f(h5_root_id,'iw',h5_w_id,ierr)
    do i=1,nepts
      write(c8,'(I8.8)')i
      call h5gcreate_f(h5_w_id,c8,h5_iw_id,ierr)
      call h5gclose_f(h5_iw_id,ierr)
    enddo
    call h5gclose_f(h5_w_id,ierr)
    call h5fclose_f(h5_root_id,ierr)
    call write_integer(nepts,1,trim(fname),'/parameters','nepts')
    call write_integer(lr_igq0,1,trim(fname),'/parameters','lr_igq0')
    call write_integer(gshme1,1,trim(fname),'/parameters','gshme1')
    call write_integer(gshme2,1,trim(fname),'/parameters','gshme2')
    call write_integer(gvecme1,1,trim(fname),'/parameters','gvecme1')
    call write_integer(gvecme2,1,trim(fname),'/parameters','gvecme2')
    call write_integer(ngvecme,1,trim(fname),'/parameters','ngvecme')
    call write_integer(spin_me,1,trim(fname),'/parameters','spin_me')
    call write_integer(nspin_chi0,1,trim(fname),'/parameters','nspin_chi0')
    call write_real8(vq0l,3,trim(fname),'/parameters','vq0l')
    call write_real8(vq0rl,3,trim(fname),'/parameters','vq0rl')
    call write_real8(vq0c,3,trim(fname),'/parameters','vq0c')
    call write_real8(vq0rc,3,trim(fname),'/parameters','vq0rc')
    call write_integer(0,1,trim(fname),'/parameters','ie1')
  else
    call read_integer(ie1,1,trim(fname),'/parameters','ie1')
  endif
endif
call i_bcast_cart(comm_cart_110,ie1,1)
ie1=ie1+1

allocate(chi0_loc(ngvecme,ngvecme,nspin_chi0))
if (lcmplx) allocate(c_chi0_loc(ngvecme,ngvecme,nspin_chi0))

if (wproc) then
  write(150,*)
  write(150,'("Starting chi0 summation")')
  write(150,'("  first energy point : ",I4)')ie1
  call flushifc(150)
endif
do ie=ie1,nepts
  call timer_reset(1)
  call timer_start(1)
  if (lcmplx) then
    c_chi0_loc=cmplx(0.d0,0.d0)
  else 
    chi0_loc=dcmplx(0.d0,0.d0)
  endif
  j=0
  do ikloc=1,nkptnr_loc
    ik=ikptnrloc(mpi_x(1),1)+ikloc-1
    call idxbos(nme(ikloc),mpi_dims(2),mpi_x(2)+1,idx0,bs)
    i1=idx0+1
    i2=idx0+bs
    j=j+bs
    do i=i1,i2
      if (nspin_chi0.eq.1) then
        ispn=1
      else
        ispn=ime(3,i,ikloc)
      endif
      wt=docc(i,ikloc)/(evalsvnr(ime(1,i,ikloc),ik) - &
        evalsvnr(ime(2,i,ikloc),idxkq(1,ikloc))+lr_w(ie))
      c_wt=cmplx(wt)
      if (lcmplx) then
        call cgerc(ngvecme,ngvecme,c_wt,c_me(1,i,ikloc),1,c_me(1,i,ikloc),1, &
          c_chi0_loc(1,1,ispn),ngvecme)
      else
        call zgerc(ngvecme,ngvecme,wt,me(1,i,ikloc),1,me(1,i,ikloc),1, &
          chi0_loc(1,1,ispn),ngvecme)
      endif
    enddo !i
  enddo !ikloc
  if (lcmplx) then
    if (mpi_dims(2).gt.1) then
      call s_reduce_cart(comm_cart_010,.false.,c_chi0_loc,2*ngvecme*ngvecme*nspin_chi0)
    endif
    if (root_cart((/0,1,0/)).and.mpi_dims(1).gt.1) then
      call s_reduce_cart(comm_cart_100,.false.,c_chi0_loc,2*ngvecme*ngvecme*nspin_chi0)
    endif
  else
    if (mpi_dims(2).gt.1) then
      call d_reduce_cart(comm_cart_010,.false.,chi0_loc,2*ngvecme*ngvecme*nspin_chi0)
    endif
    if (root_cart((/0,1,0/)).and.mpi_dims(1).gt.1) then
      call d_reduce_cart(comm_cart_100,.false.,chi0_loc,2*ngvecme*ngvecme*nspin_chi0)
    endif
  endif
  if (root_cart((/1,1,0/))) then
    if (lcmplx) chi0_loc=dcmplx(c_chi0_loc)
    chi0_loc=chi0_loc/nkptnr/omega
    write(path,'("/iw/",I8.8)')ie
    call write_real8(lr_w(ie),2,trim(fname),trim(path),'w')
    call write_real8_array(chi0_loc,4,(/2,ngvecme,ngvecme,nspin_chi0/), &
      trim(fname),trim(path),'chi0')
    call rewrite_integer(ie,1,trim(fname),'/parameters','ie1')
  endif
  call timer_stop(1)
  if (wproc) then
    write(150,'("energy point ",I4," done in ",F8.2," seconds, ",F8.2," MB/s")') &
      ie,timer(1,2),(16.0*j*ngvecme**2)/1024/1024/timer(1,2)
    call flushifc(150)
  endif
  call barrier(comm_cart_110)
enddo !ie

call barrier(comm_cart)

deallocate(lr_w)
deallocate(idxkq)
deallocate(nme)
deallocate(ime)
deallocate(docc)
deallocate(me)
deallocate(chi0_loc)
if (lcmplx) then
  deallocate(c_me)
  deallocate(c_chi0_loc)
endif

return
end
#endif
