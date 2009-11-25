#ifdef _HDF5_
subroutine response_me(ivq0m,wfsvmtloc,wfsvitloc,ngknr,igkignr,occsvnr,&
  evalsvnr,pmat)
use modmain
use hdf5
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
! q-vector in k-mesh coordinates
integer, intent(in) :: ivq0m(3)
complex(8), intent(in) :: wfsvmtloc(lmmaxvr,nrfmax,natmtot,nspinor,nstsv,*)
complex(8), intent(in) :: wfsvitloc(ngkmax,nspinor,nstsv,*)
integer, intent(in) :: ngknr(*)
integer, intent(in) :: igkignr(ngkmax,*)
real(8), intent(in) :: occsvnr(nstsv,nkptnr)
real(8), intent(in) :: evalsvnr(nstsv,nkptnr)
complex(8), intent(in) :: pmat(3,nstsv,nstsv,*)
! G-vector which brings q to first BZ
integer vgq0l(3)
! G+q vectors in Cart.coord.
real(8), allocatable :: vgq0c(:,:)
! length of G+q vectors
real(8), allocatable :: gq0(:)
! theta and phi angles of G+q vectors
real(8), allocatable :: tpgq0(:,:)
! sperical harmonics of G+q vectors
complex(8), allocatable :: ylmgq0(:,:)
! structure factor for G+q vectors
complex(8), allocatable :: sfacgq0(:,:)

! allocatable arrays
integer, allocatable :: igkignr2(:)
complex(8), allocatable :: wfsvmt2(:,:,:,:,:)
complex(8), allocatable :: wfsvit2(:,:,:)
complex(8), allocatable :: wann_c2(:,:)
integer, allocatable :: igfft1(:,:)


integer n1,n2,i1,i2,i3,itr,n,ist1,ist2
real(8) vtrc(3)

complex(8), allocatable :: gvit(:,:,:)

integer i,j,ik,jk,ig,ikstep,ierr,sz,ikloc,complete
integer ngknr2
real(8) vkq0l(3)
integer ivg1(3),ivg2(3)
real(8), allocatable :: uuj(:,:,:,:,:,:,:)
complex(4), allocatable :: gu(:,:,:)
integer, allocatable :: igu(:,:,:,:)
integer, allocatable :: ngu(:,:)
integer ngumax
complex(8) zt1

integer lmaxexp
integer lmmaxexp

character*100 :: qnm,fname
logical wproc,l1

integer, allocatable :: igishell(:)
integer, allocatable :: ishellng(:,:)
integer ngsh,gshq0

logical lwritemek
integer ik1
logical exist

! HDF5
integer(hid_t) h5_root_id
integer(hid_t) h5_kpoints_id
integer(hid_t) h5_kpoint_id
integer(hid_t) h5_tmp_id
character*8 c8

! external functions
real(8), external :: r3taxi
logical, external :: root_cart
integer, external :: iknrglob2
logical, external :: wann_diel

! comment:
! the subroutine computes <psi_{n,k}|e^{-i(G+q)x}|psi_{n',k+q}> 
!
call qname(ivq0m,qnm)
qnm="./"//trim(qnm)//"/"//trim(qnm)
wproc=.false.
if (root_cart((/1,1,0/))) then
  wproc=.true.
  fname=trim(qnm)//"_ME.OUT"
  open(150,file=trim(fname),form='formatted',status='replace')
endif

complete=0
fname=trim(qnm)//"_me.hdf5"
if (root_cart((/1,1,0/))) then
  inquire(file=trim(fname),exist=exist)
  if (exist) then
    call read_integer(complete,1,trim(fname),'/parameters','complete')
  endif
endif
call i_bcast_cart(comm_cart_110,complete,1)
if (complete.eq.1) goto 30


lmaxexp=10 !lmaxvr
lmmaxexp=(lmaxexp+1)**2

! q-vector in lattice coordinates
do i=1,3
  vq0l(i)=1.d0*(ivq0m(i))/ngridk(i)+1d-12
enddo
! find G-vector which brings q0 to first BZ
vgq0l(:)=floor(vq0l(:))
! find G-shell for a given q-vector
allocate(igishell(ngvec))
allocate(ishellng(ngvec,2))
call getgshells(ngsh,igishell,ishellng)
gshq0=igishell(ivgig(vgq0l(1),vgq0l(2),vgq0l(3)))

if (wproc) then
  write(150,*)
  write(150,'("G-shell of a given q-vector : ",I4)')gshq0
endif

if (gshq0.lt.gshme1) then
  if (wproc) then
    write(150,*)
    write(150,'("Warning: minimum number of G-shells was changed from ",&
      &I4," to ",I4)')gshme1,gshq0
  endif
  gshme1=gshq0
endif
if (gshq0.gt.gshme2) then
  if (wproc) then
    write(150,*)
    write(150,'("Warning: maximum number of G-shells was changed from ",&
      &I4," to ",I4)')gshme2,gshq0
  endif
  gshme2=gshq0
endif
! test if G-shells are closed
i=ishellng(gshme1,2)
j=ishellng(gshme2,2)
if (abs(gc(i)-gc(i+1)).lt.epslat.or.abs(gc(j)-gc(j+1)).lt.epslat) then
  write(*,*)
  write(*,'("Bug(response_me): G-shells are not closed")')
  write(*,*)
  call pstop
endif
if (gshme1.eq.1) then
  gvecme1=1
else
  gvecme1=ishellng(gshme1-1,2)+1
endif
gvecme2=ishellng(gshme2,2)
ngvecme=gvecme2-gvecme1+1

if (lscalar) then
  if (wproc) then
    write(150,*)
    write(150,'("Scalar calculation")')
  endif
  gvecme1=ivgig(vgq0l(1),vgq0l(2),vgq0l(3))
  gvecme2=gvecme1
  ngvecme=1
endif

if (wproc) then
  write(150,*)
  write(150,'("G-shell limits      : ",2I4)')gshme1,gshme2
  write(150,'("G-vector limits     : ",2I4)')gvecme1,gvecme2
  write(150,'("number of G-vectors : ",I4)')ngvecme   
  call flushifc(150)
endif
deallocate(igishell)
deallocate(ishellng)

if (wproc) then
  write(150,*)
  write(150,'("Calculation of matrix elements:")')
  write(150,'("  <n,k|e^{-i(G+q)x}|n'',k+q>")')
endif

allocate(idxkq(2,nkptnr))
allocate(vgq0c(3,ngvecme))
allocate(gq0(ngvecme))
allocate(tpgq0(2,ngvecme))
allocate(sfacgq0(ngvecme,natmtot))
allocate(ylmgq0(lmmaxexp,ngvecme))

! reduce q0 vector to first BZ
vq0rl(:)=vq0l(:)-vgq0l(:)

! check if we have enough G-shells to bring q-vector back to first BZ
do ig=1,ngvecme
  if (sum(abs(vgq0l(:)-ivg(:,ig+gvecme1-1))).eq.0) then
    lr_igq0=ig+gvecme1-1
    goto 20
  endif
enddo
write(*,*)
write(*,'("Bug(response_me): no G-vector to reduce q-vector to first BZ")')
write(*,'("  ngvecme : ",I4)')ngvecme
write(*,'("  vq0l : ",3G18.10)')vq0l
write(*,'("  vgq0l : ",3G18.10)')vgq0l
write(*,*)
call pstop
20 continue

! get Cartesian coordinates of q-vector and reduced q-vector
call r3mv(bvec,vq0l,vq0c)
call r3mv(bvec,vq0rl,vq0rc)

! write some info  
if (wproc) then
  write(150,*)
  write(150,'("q-vector (lat.coord.)                        : ",&
    & 3G18.10)')vq0l
  write(150,'("q-vector (Cart.coord.) [a.u.]                : ",&
    & 3G18.10)')vq0c
  write(150,'("q-vector length [a.u.]                       : ",&
    & G18.10)')sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)
  write(150,'("q-vector length [1/A]                        : ",&
    & G18.10)')sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)/au2ang
  write(150,'("G-vector to reduce q to first BZ (lat.coord.): ",&
    & 3I4)')vgq0l
  write(150,'("index of G-vector                            : ",&
    & I4)')lr_igq0
  write(150,'("reduced q-vector (lat.coord.)                : ",&
    & 3G18.10)')vq0rl
  write(150,'("reduced q-vector (Cart.coord.) [a.u.]        : ",&
    & 3G18.10)')vq0rc
  call flushifc(150)
endif

allocate(igfft1(ngvecme,nkptnr))
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
! search for new fft indexes
  do ig=1,ngvecme
    ivg2(:)=ivg(:,ig)+ivg1(:)
    igfft1(ig,ik)=igfft(ivgig(ivg2(1),ivg2(2),ivg2(3)))
  enddo
enddo

! hack for q=0
if (ivq0m(1).eq.0.and.ivq0m(2).eq.0.and.ivq0m(3).eq.0) then
  !vq0c=(/0.0112d0,0.0223d0,0.0331d0/)
  vq0c=(/-0.05685621179E-01,  0.05685621179E-01,  0.05685621179E-01/)
  vq0rc=vq0c
endif

! setup n,n' stuff
call timer_reset(1)
call timer_start(1)
if (spinpol) then
  allocate(spinor_ud(2,nstsv,nkptnr))
  spinor_ud=0
  do ikloc=1,nkptnr_loc
    ik=iknrglob2(ikloc,mpi_x(1))
    do j=1,nstsv
      zt1=sum(abs(evecsvloc(1:nstfv,j,ikloc)))
      if (abs(zt1).gt.1d-10) spinor_ud(1,j,ik)=1
      zt1=sum(abs(evecsvloc(nstfv+1:nstsv,j,ikloc)))
      if (abs(zt1).gt.1d-10) spinor_ud(2,j,ik)=1
    enddo
  enddo
  call i_reduce_cart(comm_cart_100,.true.,spinor_ud,2*nstsv*nkptnr)
  call i_bcast_cart(comm_cart_010,spinor_ud,2*nstsv*nkptnr)
endif
call getmeidx(.true.,occsvnr,evalsvnr,wproc)
#ifdef _MPI_
call mpi_allreduce(nmegqblhmax,i,1,MPI_INTEGER,MPI_MAX,comm_cart_100,ierr)
nmegqblhmax=i
#endif
allocate(nmegqblh(nkptnr_loc))
allocate(bmegqblh(2,nmegqblhmax,nkptnr_loc))
call getmeidx(.false.,occsvnr,evalsvnr,wproc)
call timer_stop(1)
if (wproc) then
  write(150,*)
  write(150,'("Maximum number of interband transitions: ",I5)')nmegqblhmax
  write(150,'("Done in ",F8.2," seconds")')timer(1,2)
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
    write(150,'("Number of WF trainstions : ",I4)')nmegqwan
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

! generate G+q' vectors, where q' is reduced q-vector
i=0
do ig=gvecme1,gvecme2
  i=i+1
  vgq0c(:,i)=vgc(:,ig)+vq0rc(:)
! get spherical coordinates and length of G+q'
  call sphcrd(vgq0c(:,i),gq0(i),tpgq0(:,i))
! generate spherical harmonics for G+q'
  call genylm(lmaxexp,tpgq0(:,i),ylmgq0(:,i))
enddo

! generate structure factor for G+q' vectors
call gensfacgp(ngvecme,vgq0c,ngvecme,sfacgq0)

if (root_cart((/1,1,0/))) then
  call h5fcreate_f(trim(fname),H5F_ACC_TRUNC_F,h5_root_id,ierr)
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
  call write_integer(nkptnr,1,trim(fname),'/parameters','nkptnr')
  call write_integer(nmegqblhmax,1,trim(fname),'/parameters','nmegqblhmax')
  call write_integer(lr_igq0,1,trim(fname),'/parameters','lr_igq0')
  call write_integer(gshme1,1,trim(fname),'/parameters','gshme1')
  call write_integer(gshme2,1,trim(fname),'/parameters','gshme2')
  call write_integer(gvecme1,1,trim(fname),'/parameters','gvecme1')
  call write_integer(gvecme2,1,trim(fname),'/parameters','gvecme2')
  call write_integer(ngvecme,1,trim(fname),'/parameters','ngvecme')
  call write_integer(nspinor,1,trim(fname),'/parameters','nspinor')
  call write_real8(vq0l,3,trim(fname),'/parameters','vq0l')
  call write_real8(vq0rl,3,trim(fname),'/parameters','vq0rl')
  call write_real8(vq0c,3,trim(fname),'/parameters','vq0c')
  call write_real8(vq0rc,3,trim(fname),'/parameters','vq0rc')
  complete=0
  call write_integer(complete,1,trim(fname),'/parameters','complete')
  if (wannier) then
    call write_real8(wann_occ,nwann,trim(fname),'/wannier','wann_occ')
  endif
endif
if ((.not.lsfio).and.root_cart((/0,1,0/))) then
  do ikloc=1,nkptnr_loc
    call idxglob(nkptnr,mpi_dims(1),mpi_x(1)+1,ikloc,ik)
    write(c8,'(I8.8)')ik
    fname=trim(qnm)//"_me_k_"//c8//".hdf5"
    call h5fcreate_f(trim(fname),H5F_ACC_TRUNC_F,h5_root_id,ierr)
    call h5gcreate_f(h5_root_id,'kpoints',h5_kpoints_id,ierr)
    call h5gcreate_f(h5_kpoints_id,c8,h5_kpoint_id,ierr)
    call h5gclose_f(h5_kpoint_id,ierr)
    call h5gclose_f(h5_kpoints_id,ierr)
    call h5fclose_f(h5_root_id,ierr)
  enddo
endif

call timer_reset(1)
call timer_start(1)
if (wproc) then
  write(150,*)
  write(150,'("Calculating radial integrals")')
  write(150,'("  maximum number of radial functions : ",I4)')nrfmax
endif
allocate(uuj(0:lmaxvr,0:lmaxvr,0:lmaxexp,nrfmax,nrfmax,natmtot,ngvecme))
call calc_uuj(uuj,lmaxexp,gq0)
call timer_stop(1)
if (wproc) then
  write(150,'("Done in ",F8.2," seconds")')timer(1,2)
  call flushifc(150)
endif

call timer_reset(1)
call timer_start(1)
call getgu(.true.,lmaxexp,uuj,ylmgq0,sfacgq0,ngumax,ngu,gu,igu)
allocate(ngu(natmtot,ngvecme))
allocate(gu(ngumax,natmtot,ngvecme))
allocate(igu(4,ngumax,natmtot,ngvecme))
call getgu(.false.,lmaxexp,uuj,ylmgq0,sfacgq0,ngumax,ngu,gu,igu)
call timer_stop(1)
sz=natmtot*ngvecme+ngumax*natmtot*ngvecme+4*ngumax*natmtot*ngvecme
sz=8*sz/1024/1024
if (wproc) then
  write(150,*)
  write(150,'("Maximum number of Gaunt-like coefficients : ",I8)')ngumax
  write(150,'("  size (MB) : ",I6)')sz
  write(150,'("Done in ",F8.2," seconds")')timer(1,2)
  call flushifc(150)
endif
deallocate(uuj)

allocate(gvit(intgv(1,1):intgv(1,2),intgv(2,1):intgv(2,2),intgv(3,1):intgv(3,2)))
call gengvit(gvit)

sz=16*ngvecme*nmegqblhmax*nkptnrloc(0)/1024/1024
if (sz.lt.150) then
  if (wproc) then
    write(150,*)
    write(150,'("Size of matrix elements (Mb): ",I6)')sz
    call flushifc(150)
  endif
  lwritemek=.false.
else
  lwritemek=.true.
endif
if (wannier) lwritemek=.true.
if (lwritemek) then
  allocate(megqblh(ngvecme,nmegqblhmax,1))
else
  allocate(megqblh(ngvecme,nmegqblhmax,nkptnr_loc))
endif
allocate(wfsvmt2(lmmaxvr,nrfmax,natmtot,nspinor,nstsv))
allocate(wfsvit2(ngkmax,nspinor,nstsv))
allocate(igkignr2(ngkmax))
if (wannier) then
  allocate(wann_c2(nwann,nstsv))
endif

if (wproc) then
  write(150,*)
  write(150,'("Starting k-point loop")')
  call flushifc(150)
endif
do ikstep=1,nkptnrloc(0)
  call timer_reset(1)
  call timer_reset(2)
  call timer_reset(3)
  call timer_reset(4)
  call timer_reset(5)
  call timer_reset(6)
  if (wproc) then
    write(150,'("k-step ",I4," out of ",I4)')ikstep,nkptnrloc(0)
    call flushifc(150)
  endif
! transmit wave-functions
  call timer_start(1)
  call wfkq(ikstep,wfsvmtloc,wfsvitloc,ngknr,igkignr,wfsvmt2, &
    wfsvit2,ngknr2,igkignr2,wann_c2)
  call barrier(comm_cart_110)
  if (wproc) then
    write(150,'("  wave-functions are distributed")')
    call flushifc(150)
  endif
  call timer_stop(1)
! compute matrix elements  
  call timer_start(2)
  if (ikstep.le.nkptnrloc(mpi_x(1))) then
    call idxglob(nkptnr,mpi_dims(1),mpi_x(1)+1,ikstep,ik)
    if (lwritemek) then
      ik1=1
    else
      ik1=ikstep
    endif
    megqblh(:,:,ik1)=zzero
    if (.not.lmeoff) then
! calculate muffin-tin contribution for all combinations of n,n'
      call timer_start(4)
      call megqblhmt(nmegqblh(ikstep),bmegqblh(1,1,ikstep),               &
        wfsvmtloc(1,1,1,1,1,ikstep),wfsvmt2,ngumax,ngu,gu,igu, &
        megqblh(1,1,ik1),ik,idxkq(1,ik))
      call timer_stop(4)
! calculate interstitial contribution for all combinations of n,n'
      call timer_start(5)
      call megqblhit(nmegqblh(ikstep),bmegqblh(1,1,ikstep),ngknr(ikstep), &
        ngknr2,igkignr(1,ikstep),igkignr2,idxkq(2,ik),         &
        wfsvitloc(1,1,1,ikstep),wfsvit2,megqblh(1,1,ik1),gvit,ik,idxkq(1,ik))
      call timer_stop(5)
    else
      megqblh(1:ngvecme,1:nmegqblh(ikstep),ik1)=zone
    endif
! hack for q=0
    if (ivq0m(1).eq.0.and.ivq0m(2).eq.0.and.ivq0m(3).eq.0) then
      megqblh(1,:,ik1)=zzero
      do i=1,nmegqblh(ikstep)
        ist1=bmegqblh(1,i,ikstep)
        ist2=bmegqblh(2,i,ikstep)
        if (ist1.eq.ist2) megqblh(1,i,ik1)=zone
        megqblh(1,i,ik1)=megqblh(1,i,ik1)-&
          dot_product(vq0rc(:),pmat(:,ist1,ist2,ikstep))/&
          (evalsvnr(ist1,ik)-evalsvnr(ist2,ik)+swidth)          
      enddo
    endif
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
                  wann_c2(n2,ist2)*megqblh(ig,i,ik1)*zt1
            enddo !i
          enddo !n
        enddo !ig      
      enddo !itr
    endif !wannier
  endif ! ikstep.le.nkptnrloc(mpi_x(1))
  call timer_stop(2)
! write matrix elements
  if (lwritemek) then
    call timer_start(3)
! only the plane of (k,q) processors will write
    if (root_cart((/0,1,0/))) then
! if writing to one file
      if (lsfio) then
        do i=0,mpi_dims(1)-1
          if (mpi_x(1).eq.i.and.ikstep.le.nkptnr_loc) then
            call writeme(ikstep,fname,megqblh,wann_c2,pmat(1,1,1,ikstep))
          endif
          call barrier(comm_cart_100)
        enddo   
      else
        if (ikstep.le.nkptnr_loc) then
          write(fname,'("_me_k_",I8.8)')ik
          fname=trim(qnm)//trim(fname)//".hdf5"
          call writeme(ikstep,fname,megqblh,wann_c2,pmat(1,1,1,ikstep))
        endif
      endif !lsfio
    endif !root_cart
    call timer_stop(3)
  endif !lwritemek
  if (wproc) then
    write(150,'("  time (seconds)")')
    write(150,'("    send/recv      : ",F8.2)')timer(1,2)
    write(150,'("    matrix elements")')
    write(150,'("      muffin-tins  (Bloch basis) : ",F8.2)')timer(4,2)
    write(150,'("      interstitial (Bloch basis) : ",F8.2)')timer(5,2)
    write(150,'("      writing      : ",F8.2)')timer(3,2)
    write(150,'("    total for kpt  : ",F8.2)')timer(2,2)+timer(3,2)
    write(150,'("  speed (me/sec)   : ",F10.2)') &
      ngvecme*nmegqblh(ikstep)/(timer(2,2)+timer(3,2))
    call flushifc(150)
  endif
enddo !ikstep

if (wannier) then
! sum over all k-points to get <n,T=0|e^{-i(G+q)x|n',T'>
  if (root_cart((/0,1,0/)).and.mpi_dims(1).gt.1) then
    call d_reduce_cart(comm_cart_100,.false.,megqwan,2*nmegqwan*ntrmegqwan*ngvecme)
  endif
  megqwan=megqwan/nkptnr
  if (root_cart((/1,1,0/))) then
    fname=trim(qnm)//"_me.hdf5"
    call write_integer(ntrmegqwan,1,trim(fname),'/wannier','ntrmegqwan')
    call write_integer(nmegqwan,1,trim(fname),'/wannier','nmegqwan')
    call write_real8_array(megqwan,4,(/2,nmegqwan,ntrmegqwan,ngvecme/), &
      trim(fname),'/wannier','megqwan')
    call write_integer_array(itrmegqwan,2,(/3,ntrmegqwan/),trim(fname),'/wannier','itrmegqwan')
    call write_integer_array(bmegqwan,2,(/2,nmegqwan/),trim(fname),'/wannier','bmegqwan')
!    open(170,file=trim(trim(qnm)//"_megqwan.txt"),status='replace',form='formatted')
!    write(170,*)sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2),&
!      abs(1.d0-megqwan(1,1,1))/sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2),abs(megqwan(1,1,1))
!    close(170)
  endif
endif

if (.not.lwritemek) then
  if (wproc) then
    write(150,*)
    write(150,'("Writing matrix elements")')
  endif
  call timer_start(3)
  if (root_cart((/0,1,0/))) then
    if (lsfio) then
      do i=0,mpi_dims(1)-1
        if (mpi_x(1).eq.i) then
          do ikstep=1,nkptnr_loc
            call writeme(ikstep,fname,megqblh(1,1,ikstep),wann_c2,pmat(1,1,1,ikstep))
          enddo
        endif
        call barrier(comm_cart_100)
      enddo !i
    else
      do ikstep=1,nkptnr_loc
        call idxglob(nkptnr,mpi_dims(1),mpi_x(1)+1,ikstep,ik)
        write(fname,'("_me_k_",I8.8)')ik
        fname=trim(qnm)//trim(fname)//".hdf5"
        call writeme(ikstep,fname,megqblh(1,1,ikstep),wann_c2,pmat(1,1,1,ikstep))
      enddo
    endif !lsfio
  endif !root_cart
  call timer_stop(3)
  if (wproc) write(150,'(" Done in : ",F8.2)')timer(3,2)
endif

deallocate(gvit)
deallocate(wfsvmt2)
deallocate(wfsvit2)
deallocate(igkignr2)
deallocate(megqblh)
deallocate(nmegqblh)
deallocate(bmegqblh)
deallocate(idxkq)
deallocate(vgq0c)
deallocate(gq0)
deallocate(tpgq0)
deallocate(sfacgq0)
deallocate(ylmgq0) 
deallocate(ngu)
deallocate(gu)
deallocate(igu)
if (spinpol) then
  deallocate(spinor_ud)
endif
if (wannier) then
  deallocate(wann_c2)
  deallocate(bmegqwan)
  deallocate(itrmegqwan)
  deallocate(megqwan)
endif
deallocate(igfft1)

if (root_cart((/1,1,0/))) then
  complete=1
  call rewrite_integer(complete,1,trim(fname),'/parameters','complete')
endif

call barrier(comm_cart_110)

30 continue
if (wproc) then
  write(150,*)
  write(150,'("Done.")')
  call flushifc(150)
endif

return
end
#endif
