#ifdef _HDF5_
subroutine genmegq(ivq0m,wfsvmtloc,wfsvitloc,ngknr,igkignr,pmat)
use modmain
use hdf5
implicit none
! arguments
! q-vector in k-mesh coordinates
integer, intent(in) :: ivq0m(3)
complex(8), intent(in) :: wfsvmtloc(lmmaxvr,nrfmax,natmtot,nspinor,nstsv,nkptnrloc)
complex(8), intent(in) :: wfsvitloc(ngkmax,nspinor,nstsv,nkptnrloc)
integer, intent(in) :: ngknr(nkptnrloc)
integer, intent(in) :: igkignr(ngkmax,nkptnrloc)
complex(8), intent(in) :: pmat(3,nstsv,nstsv,nkptnrloc)

! allocatable arrays
integer, allocatable :: igkignr2(:)
complex(8), allocatable :: wfsvmt2(:,:,:,:,:)
complex(8), allocatable :: wfsvit2(:,:,:)
complex(8), allocatable :: wann_c2(:,:)

integer n1,n2,i1,i2,i3,itr,n,ist1,ist2,ig1,ig2
real(8) vtrc(3)

integer i,j,ik,jk,ig,ikstep,ierr,sz,ikloc,complete
integer ngknr2
real(8) vkq0l(3)
integer ivg1(3),ivg2(3)
complex(8), allocatable :: gntuju(:,:,:)
integer, allocatable :: igntuju(:,:,:,:)
integer, allocatable :: ngntuju(:,:)
integer ngntujumax
complex(8) zt1
integer nkstep

complex(8), allocatable :: chi0w(:,:)

integer lmaxexp
integer lmmaxexp

character*100 :: qnm,fout,fme,fmek
logical l1

integer, allocatable :: igishell(:)
integer, allocatable :: ishellng(:,:)
integer ngsh,gshq0

logical exist

complex(8), allocatable :: krnl(:,:)
complex(8), allocatable :: krnl_scr(:,:)
real(8), allocatable :: vcgq(:)
complex(8), allocatable :: epsilon(:,:)
complex(8), allocatable :: uscrn(:,:)
complex(8), allocatable :: ubare(:,:)



! HDF5
integer(hid_t) h5_root_id
integer(hid_t) h5_kpoints_id
integer(hid_t) h5_kpoint_id
integer(hid_t) h5_tmp_id
character*8 c8

! external functions
real(8), external :: r3taxi
logical, external :: wann_diel

! comment:
! the subroutine computes <psi_{n,k}|e^{-i(G+q)x}|psi_{n',k+q}> 
!
! quck turnaround for cRPA: compute chi0 and chi here; no writing 
! 
! switch write_megq_file controls the reading and writing of ME file
! when we write ME we have two choices: write to single file or write
!  to multiple files

! maximum l for exponent expansion
lmaxexp=lmaxvr+2
lmmaxexp=(lmaxexp+1)**2

call qname(ivq0m,qnm)
qnm="./"//trim(qnm)//"/"//trim(qnm)
wproc=.false.
if (mpi_grid_root((/dim_k,dim2/))) then
  wproc=.true.
  fout=trim(qnm)//"_ME.OUT"
  open(150,file=trim(fout),form='formatted',status='replace')
endif

complete=0
fme=trim(qnm)//"_me.hdf5"
if (mpi_grid_root((/dim_k,dim2/))) then
  inquire(file=trim(fme),exist=exist)
  if (exist) then
    call read_integer(complete,1,trim(fme),'/parameters','complete')
  endif
endif
call mpi_grid_bcast(complete,dims=(/dim_k,dim2/))
if (complete.eq.1) goto 30

if (wproc) then
  write(150,*)
  write(150,'("Calculation of matrix elements:")')
  write(150,'("  <n,k|e^{-i(G+q)x}|n'',k+q>")')
endif

call init_g_q_gq(ivq0m,lmaxexp,lmmaxexp)

allocate(idxkq(2,nkptnr))
! find k+q and reduce them to first BZ (this is required to utilize the 
!   periodical property of Bloch-states: |k>=|k+K>, where K is any vector 
!   of the reciprocal lattice)
do ik=1,nkptnr
! k+q vector
  vkq0l(:)=vklnr(:,ik)+vq0rl(:)+1d-12
! K vector
  ivg1(:)=floor(vkq0l(:))
! reduced k+q vector: k'=k+q-K
  vkq0l(:)=vkq0l(:)-ivg1(:)
! search for index of reduced k+q vector 
  do jk=1,nkptnr
    if (r3taxi(vklnr(:,jk),vkq0l).lt.epslat) then
      idxkq(1,ik)=jk
      goto 10
    endif
  enddo
  write(*,*)
  write(*,'("Error(response_me): index of reduced k+q point is not found")')
  write(*,'(" index of k-point: ",I4)')ik
  write(*,'(" K-vector: ",3I4)')ivg1
  write(*,'(" reduced k+q vector: ",3G18.10)')vkq0l
  write(*,'(" check original q-vector coordinates")')
  write(*,*)
  call pstop
10 continue
  idxkq(2,ik)=ivgig(ivg1(1),ivg1(2),ivg1(3))
enddo

! hack for q=0
!if (ivq0m(1).eq.0.and.ivq0m(2).eq.0.and.ivq0m(3).eq.0) then
!  !vq0c=(/0.0112d0,0.0223d0,0.0331d0/)
!  vq0c=(/-0.05685621179E-01,  0.05685621179E-01,  0.05685621179E-01/)
!  vq0rc=vq0c
!endif

! setup n,n' stuff
call timer_start(1,reset=.true.)
if (spinpol) then
  allocate(spinor_ud(2,nstsv,nkptnr))
  spinor_ud=0
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do j=1,nstsv
      zt1=sum(abs(evecsvloc(1:nstfv,j,ikloc)))
      if (abs(zt1).gt.1d-10) spinor_ud(1,j,ik)=1
      zt1=sum(abs(evecsvloc(nstfv+1:nstsv,j,ikloc)))
      if (abs(zt1).gt.1d-10) spinor_ud(2,j,ik)=1
    enddo
  enddo
  call mpi_grid_reduce(spinor_ud(1,1,1),2*nstsv*nkptnr,dims=(/dim_k/),all=.true.)
endif
call getmeidx(.true.)
call mpi_grid_reduce(nmegqblhmax,dims=(/dim_k/),all=.true.,op=op_max)
allocate(nmegqblh(nkptnrloc))
allocate(bmegqblh(2,nmegqblhmax,nkptnrloc))
call getmeidx(.false.)
call timer_stop(1)
if (wproc) then
  write(150,*)
  write(150,'("Maximum number of interband transitions: ",I5)')nmegqblhmax
  write(150,'("Done in ",F8.2," seconds")')timer_get_value(1)
  call flushifc(150)
endif

if (wannier) then
  allocate(bmegqwan(2,nwann*nwann))
  nmegqwan=0
  do n1=1,nwann
    do n2=1,nwann
      l1=.false.
! for integer occupancy numbers take only transitions between occupied and empty bands
      if (wann_diel().and.(abs(wann_occ(n1)-wann_occ(n2)).gt.1d-8)) l1=.true.
! for fractional occupancies or cRPA calculation take all transitions
      if (.not.wann_diel().or.crpa) l1=.true.
      if (l1) then
        nmegqwan=nmegqwan+1
        bmegqwan(1,nmegqwan)=n1
        bmegqwan(2,nmegqwan)=n2
      endif
    enddo
  enddo
  if (wproc) then
    write(150,*)
    write(150,'("Number of WF transitions : ",I4)')nmegqwan
  endif
! list of translations
  ntrmegqwan=(2*lr_maxtr+1)**3
  allocate(itrmegqwan(3,ntrmegqwan))
  i=0
  do i1=-lr_maxtr,lr_maxtr
    do i2=-lr_maxtr,lr_maxtr
      do i3=-lr_maxtr,lr_maxtr
        i=i+1
        itrmegqwan(:,i)=(/i1,i2,i3/)
      enddo
    enddo
  enddo
  allocate(megqwan(nmegqwan,ntrmegqwan,ngvecme))
  megqwan=zzero
endif

if (mpi_grid_root(dims=(/dim_k,dim2/)).and.write_megq_file) then
  call h5fcreate_f(trim(fme),H5F_ACC_TRUNC_F,h5_root_id,ierr)
  call h5gcreate_f(h5_root_id,'parameters',h5_tmp_id,ierr)
  call h5gclose_f(h5_tmp_id,ierr)
  if (wannier) then
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
  if (wannier) then
    call write_real8(wann_occ,nwann,trim(fme),'/wannier','wann_occ')
  endif
endif
if (mpi_grid_root(dims=(/dim_k,dim2/)).and.write_megq_file.and.split_megq_file) then
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

call getmaxgnt(lmaxexp,ngntujumax)
ngvecmeloc=mpi_grid_map(ngvecme,dim_g)

call timer_start(1,reset=.true.)
allocate(ngntuju(natmtot,ngvecmeloc))
allocate(igntuju(4,ngntujumax,natmtot,ngvecmeloc))
allocate(gntuju(ngntujumax,natmtot,ngvecmeloc))
call gengntuju(lmaxexp,ngntujumax,ngntuju,igntuju,gntuju)
call timer_stop(1)
sz=4.d0*natmtot*ngvecmeloc*(1+8*ngntujumax)/1024/1024
if (wproc) then
  write(150,*)
  write(150,'("Maximum number of Gaunt-like coefficients : ",I8)')ngntujumax
  write(150,'("  size (MB) : ",I6)')sz
  write(150,'("Done in ",F8.2," seconds")')timer_get_value(1)
  call flushifc(150)
endif

sz=16.d0*ngvecme*nmegqblhmax*nkptnrloc/1024/1024
if (wproc) then
  write(150,*)
  write(150,'("Size of matrix elements (MB): ",I6)')sz
  call flushifc(150)
endif
allocate(megqblh(ngvecme,nmegqblhmax,nkptnrloc))
allocate(wfsvmt2(lmmaxvr,nrfmax,natmtot,nspinor,nstsv))
allocate(wfsvit2(ngkmax,nspinor,nstsv))
allocate(igkignr2(ngkmax))
if (wannier) then
  allocate(wann_c2(nwann,nstsv))
endif

if (crpa) then
  allocate(chi0w(ngvecme,ngvecme))
  chi0w=zzero
endif

i=0
nkstep=mpi_grid_map(nkptnr,dim_k,x=i)
if (wproc) then
  write(150,*)
  write(150,'("Starting k-point loop")')
  call flushifc(150)
endif
do ikstep=1,nkstep
  if (wproc) then
    write(150,'("k-step ",I4," out of ",I4)')ikstep,nkstep
    call flushifc(150)
  endif
! transmit wave-functions
  call timer_start(1,reset=.true.)
  call getwfkq(ikstep,wfsvmtloc,wfsvitloc,ngknr,igkignr,wfsvmt2, &
    wfsvit2,ngknr2,igkignr2,wann_c2)
  if (wproc) then
    write(150,'("  wave-functions are distributed")')
    call flushifc(150)
  endif
  call timer_stop(1)
! compute matrix elements  
  call timer_start(2,reset=.true.)
  if (ikstep.le.nkptnrloc) then
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikstep)
    jk=idxkq(1,ik)
    megqblh(:,:,ikstep)=zzero
! calculate muffin-tin contribution for all combinations of n,n'
    call timer_start(4,reset=.true.)
    call megqblhmt(nmegqblh(ikstep),bmegqblh(1,1,ikstep),&
      wfsvmtloc(1,1,1,1,1,ikstep),wfsvmt2,ngntujumax,ngntuju,igntuju,gntuju,&
      ik,jk,megqblh(1,1,ikstep))
    call timer_stop(4)
! calculate interstitial contribution for all combinations of n,n'
    call timer_start(5,reset=.true.)
    call megqblhit(nmegqblh(ikstep),bmegqblh(1,1,ikstep),ngknr(ikstep), &
      ngknr2,igkignr(1,ikstep),igkignr2,idxkq(2,ik),         &
      wfsvitloc(1,1,1,ikstep),wfsvit2,ik,jk,megqblh(1,1,ikstep))
    call timer_stop(5)
! hack for q=0
!    if (ivq0m(1).eq.0.and.ivq0m(2).eq.0.and.ivq0m(3).eq.0) then
!      megqblh(1,:,ik1)=zzero
!      do i=1,nmegqblh(ikstep)
!        ist1=bmegqblh(1,i,ikstep)
!        ist2=bmegqblh(2,i,ikstep)
!        if (ist1.eq.ist2) megqblh(1,i,ik1)=zone
!        megqblh(1,i,ik1)=megqblh(1,i,ik1)-&
!          dot_product(vq0rc(:),pmat(:,ist1,ist2,ikstep))/&
!          (evalsvnr(ist1,ik)-evalsvnr(ist2,ik)+swidth)          
!      enddo
!    endif
! add contribution from k-point to the matrix elements of e^{-i(G+q)x} in 
!  the basis of Wannier functions
    if (wannier) then
      do itr=1,ntrmegqwan
        vtrc(:)=avec(:,1)*itrmegqwan(1,itr)+&
                avec(:,2)*itrmegqwan(2,itr)+&
                avec(:,3)*itrmegqwan(3,itr)
        zt1=exp(dcmplx(0.d0,-dot_product(vkcnr(:,ik)+vq0rc(:),vtrc(:))))
        do ig=1,ngvecme
          do n=1,nmegqwan
            n1=bmegqwan(1,n)
            n2=bmegqwan(2,n)
            do i=1,nmegqblh(ikstep)
              ist1=bmegqblh(1,i,ikstep)
              ist2=bmegqblh(2,i,ikstep)
              megqwan(n,itr,ig)=megqwan(n,itr,ig)+dconjg(wann_c(n1,ist1,ikstep))*&
                  wann_c2(n2,ist2)*megqblh(ig,i,ikstep)*zt1
            enddo !i
          enddo !n
        enddo !ig      
      enddo !itr
    endif !wannier
!    if (crpa) then
! for each k-point : sum over interband transitions
!      call sum_chi0(ikstep,ikstep,ik,idxkq(1,ik),nmegqblh(ikstep),1,nmegqblh(ikstep),&
!        lr_evalsvnr,lr_occsvnr,zi*lr_eta/ha2ev,chi0w)
!    endif
  endif ! ikstep.le.nkptnrloc(mpi_x(1))
  call timer_stop(2)
  if (wproc) then
    write(150,'("  time (seconds)")')
    write(150,'("    send/recv      : ",F8.2)')timer_get_value(1)
    write(150,'("    matrix elements")')
    write(150,'("      muffin-tins  (Bloch basis) : ",F8.2)')timer_get_value(4)
    write(150,'("      interstitial (Bloch basis) : ",F8.2)')timer_get_value(5)
    write(150,'("    total for kpt  : ",F8.2)')timer_get_value(2)
    write(150,'("  speed (me/sec)   : ",F10.2)')&
      ngvecme*nmegqblh(ikstep)/timer_get_value(2)
    call flushifc(150)
  endif
enddo !ikstep

!if (wannier) then
!! sum over all k-points to get <n,T=0|e^{-i(G+q)x|n',T'>
!  if (root_cart((/0,1,0/)).and.mpi_dims(1).gt.1) then
!    call d_reduce_cart(comm_cart_100,.false.,megqwan,2*nmegqwan*ntrmegqwan*ngvecme)
!  endif
!  megqwan=megqwan/nkptnr
!  if (root_cart((/1,1,0/)).and.lwriteme) then
!    fname=trim(qnm)//"_me.hdf5"
!    call write_integer(ntrmegqwan,1,trim(fname),'/wannier','ntrmegqwan')
!    call write_integer(nmegqwan,1,trim(fname),'/wannier','nmegqwan')
!    call write_real8_array(megqwan,4,(/2,nmegqwan,ntrmegqwan,ngvecme/), &
!      trim(fname),'/wannier','megqwan')
!    call write_integer_array(itrmegqwan,2,(/3,ntrmegqwan/),trim(fname),'/wannier','itrmegqwan')
!    call write_integer_array(bmegqwan,2,(/2,nmegqwan/),trim(fname),'/wannier','bmegqwan')
!  endif
!endif

if (write_megq_file) then
  if (wproc) then
    write(150,*)
    write(150,'("Writing matrix elements")')
  endif
  call timer_start(3,reset=.true.)
  if (mpi_grid_root(dims=(/dim2/))) then
    if (.not.split_megq_file) then
      do i=0,mpi_grid_size(1)-1
        if (mpi_grid_x(1).eq.i) then
          do ikloc=1,nkptnrloc
            call writeme(ikloc,fme,megqblh(1,1,ikloc),wann_c2,pmat(1,1,1,ikloc))
          enddo
        endif
        call mpi_grid_barrier(dims=(/dim_k/))
      enddo !i
    else
      do ikloc=1,nkptnrloc
        ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
        write(fmek,'("_me_k_",I8.8)')ik
        fmek=trim(qnm)//trim(fmek)//".hdf5"
        call writeme(ikloc,fmek,megqblh(1,1,ikloc),wann_c2,pmat(1,1,1,ikloc))
      enddo
    endif !.not.spit_megq_file
  endif !mpi_grid_root
  call timer_stop(3)
  if (wproc) write(150,'(" Done in : ",F8.2)')timer_get_value(3)
endif

!if (crpa) then
!! sum chi0 over k-points
!  if (root_cart((/0,1,0/)).and.mpi_dims(1).gt.1) then
!    call d_reduce_cart(comm_cart_100,.false.,chi0w,2*ngvecme*ngvecme)
!  endif
!  if (root_cart((/1,1,0/))) then
!	allocate(krnl(ngvecme,ngvecme))
!	allocate(krnl_scr(ngvecme,ngvecme))
!	allocate(vcgq(ngvecme))
!	allocate(epsilon(ngvecme,ngvecme))
!	krnl=zzero
!	krnl_scr=zzero
!	vcgq=0.d0
!	do ig=1,ngvecme
!	  vcgq(ig)=2*sqrt(pi)/gq0(ig)
!	  krnl(ig,ig)=vcgq(ig)**2
!	enddo !ig
!    chi0w=chi0w/nkptnr/omega
!	do ig1=1,ngvecme
!	  do ig2=1,ngvecme
!		epsilon(ig1,ig2)=-vcgq(ig1)*chi0w(ig1,ig2)*vcgq(ig2)
!	  enddo
!	  epsilon(ig1,ig1)=dcmplx(1.d0,0.d0)+epsilon(ig1,ig1)
!	enddo
!	call invzge(epsilon,ngvecme)
!	do ig1=1,ngvecme
!	  do ig2=1,ngvecme
!		krnl_scr(ig1,ig2)=vcgq(ig1)*epsilon(ig1,ig2)*vcgq(ig2)
!	  enddo
!	enddo
!  
!	allocate(imegqwan(nwann,nwann))
!	imegqwan=-1
!	do i=1,nmegqwan
!	  n1=bmegqwan(1,i)
!	  n2=bmegqwan(2,i)
!	  imegqwan(n1,n2)=i
!	enddo
!	
!	allocate(uscrn(nwann,nwann))
!	allocate(ubare(nwann,nwann))
!	uscrn=zzero
!	ubare=zzero
!	do i1=1,ngvecme
!	  do i2=1,ngvecme
!		do n1=1,nwann
!		  do n2=1,nwann
!			uscrn(n1,n2)=uscrn(n1,n2)+dconjg(megqwan(imegqwan(n1,n1),1,i1))*&
!			  krnl_scr(i1,i2)*megqwan(imegqwan(n2,n2),1,i2)
!			ubare(n1,n2)=ubare(n1,n2)+dconjg(megqwan(imegqwan(n1,n1),1,i1))*&
!			  krnl(i1,i2)*megqwan(imegqwan(n2,n2),1,i2)
!		  enddo
!		enddo
!	  enddo
!	enddo
!	uscrn=ha2ev*uscrn/omega
!	ubare=ha2ev*ubare/omega
!	fname=trim(qnm)//"_U"
!	open(170,file=trim(fname),status='replace',form='unformatted')
!	write(170)uscrn,ubare
!	close(170)
!	fname=trim(qnm)//"_U.txt"
!	open(170,file=trim(fname),status='replace',form='formatted')
!	write(170,'("Screened U matrix")')
!	write(170,'("real part")')
!	do i=1,nwann
!	  write(170,'(100F12.6)')(dreal(uscrn(i,j)),j=1,nwann)
!	enddo
!	write(170,'("imag part")')
!	do i=1,nwann
!	  write(170,'(100F12.6)')(dimag(uscrn(i,j)),j=1,nwann)
!	enddo
!	write(170,*)      
!	write(170,'("Bare U matrix")')
!	write(170,'("real part")')
!	do i=1,nwann
!	  write(170,'(100F12.6)')(dreal(ubare(i,j)),j=1,nwann)
!	enddo
!	write(170,'("imag part")')
!	do i=1,nwann
!	  write(170,'(100F12.6)')(dimag(ubare(i,j)),j=1,nwann)
!	enddo  
!	close(170)
!
!	if (ngvecme.gt.10) then
!	  n1=10
!	else
!	  n1=ngvecme
!	endif
!	fname=trim(qnm)//"_krnl.txt"
!	open(170,file=trim(fname),status='replace',form='formatted')
!	write(170,'("Screened W matrix")')
!	write(170,'("real part")')
!	do i=1,n1
!	  write(170,'(100F12.6)')(dreal(krnl_scr(i,j)),j=1,n1)
!	enddo
!	write(170,'("imag part")')
!	do i=1,n1
!	  write(170,'(100F12.6)')(dimag(krnl_scr(i,j)),j=1,n1)
!	enddo
!	close(170)
!	
!	deallocate(uscrn,ubare)
!	deallocate(krnl)
!	deallocate(krnl_scr)
!	deallocate(vcgq)
!	deallocate(epsilon)
!	deallocate(imegqwan)
!  endif
!endif !crpa

! deallocate arrays if we saved the ME file
if (write_megq_file) then
  deallocate(megqblh)
  deallocate(nmegqblh)
  deallocate(bmegqblh)
  deallocate(idxkq)
endif

deallocate(wfsvmt2)
deallocate(wfsvit2)
deallocate(igkignr2)
deallocate(ngntuju)
deallocate(gntuju)
deallocate(igntuju)
if (spinpol) then
  deallocate(spinor_ud)
endif
if (wannier) then
  deallocate(wann_c2)
  deallocate(bmegqwan)
  deallocate(itrmegqwan)
  deallocate(megqwan)
endif
if (crpa) then
  deallocate(chi0w)
endif

!if (root_cart((/1,1,0/)).and.lwriteme) then
!  complete=1
!  fname=trim(qnm)//"_me.hdf5"
!  call rewrite_integer(complete,1,trim(fname),'/parameters','complete')
!endif

!call barrier(comm_cart_110)

30 continue
if (wproc) then
  write(150,*)
  write(150,'("Done.")')
  call flushifc(150)
endif

return
end
#endif
