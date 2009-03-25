#ifdef _HDF5_
subroutine response_me(ivq0m,wfsvmtloc,wfsvitloc,ngknr,igkignr,occsvnr)
use modmain
use hdf5
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
! q-vector in k-mesh coordinates
integer, intent(in) :: ivq0m(3)
complex(8), intent(in) :: wfsvmtloc(lmmaxvr,nrfmax,natmtot,nstsv,nspinor,*)
complex(8), intent(in) :: wfsvitloc(ngkmax,nstsv,nspinor,*)
integer, intent(in) :: ngknr(*)
integer, intent(in) :: igkignr(ngkmax,*)
real(8), intent(in) :: occsvnr(nstsv,nkptnr)
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
! hack to switch off matrix elements
logical, parameter :: meoff=.false.
! indexes for fft-transform of u_{nk}^{*}u_{n'k'}exp{-iKx}, where k'=k+q-K
integer, allocatable :: igfft1(:,:)

! allocatable arrays
integer, allocatable :: igkignr2(:)
complex(8), allocatable :: wfsvmt2(:,:,:,:,:)
complex(8), allocatable :: wfsvit2(:,:,:)
complex(8), allocatable :: evecfv2(:,:,:)
complex(8), allocatable :: evecsv2(:,:)
complex(8), allocatable :: me(:,:)

integer i,j,ik,jk,ig,ikstep,ierr,sz,ikloc
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

character*100 :: qnm,fname,path

integer, allocatable :: igishell(:)
integer, allocatable :: ishellng(:,:)
integer ngsh,gshq0

! HDF5
integer(hid_t) h5_root_id
integer(hid_t) h5_kpoints_id
integer(hid_t) h5_kpoint_id
integer(hid_t) h5_tmp_id
character*8 c8

! external functions
real(8), external :: r3taxi
complex(8), external :: zfint
logical, external :: root_cart
integer, external :: iknrglob2

lmaxexp=lmaxvr
lmmaxexp=(lmaxexp+1)**2

! q-vector in lattice coordinates
do i=1,3
  vq0l(i)=1.d0*ivq0m(i)/ngridk(i)+1d-12
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
  if (.not.spinpol) write(150,'("  <n,k|e^{-i(G+q)x}|n'',k+q>")')
  if (spinpol.and.lrtype.eq.0) then
    if (spin_me.eq.1.or.spin_me.eq.3) write(150,'("  <n,k,up|e^{-i(G+q)x}|n'',k+q,up>")')
    if (spin_me.eq.2.or.spin_me.eq.3) write(150,'("  <n,k,dn|e^{-i(G+q)x}|n'',k+q,dn>")')
  endif
  if (spinpol.and.lrtype.eq.1) then
    if (spin_me.eq.1.or.spin_me.eq.3) write(150,'("  <n,k,up|e^{-i(G+q)x}|n'',k+q,dn>")')
    if (spin_me.eq.2.or.spin_me.eq.3) write(150,'("  <n,k,dn|e^{-i(G+q)x}|n'',k+q,up>")')
  endif
endif

allocate(idxkq(2,nkptnr))
allocate(vgq0c(3,ngvecme))
allocate(gq0(ngvecme))
allocate(tpgq0(2,ngvecme))
allocate(sfacgq0(ngvecme,natmtot))
allocate(ylmgq0(lmmaxexp,ngvecme))
allocate(igfft1(ngvecme,nkptnr))

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
endif
call getmeidx(.true.,occsvnr)
#ifdef _MPI_
call mpi_allreduce(nmemax,i,1,MPI_INTEGER,MPI_MAX,comm_cart_100,ierr)
nmemax=i
#endif
allocate(nme(nkptnr_loc))
allocate(ime(3,nmemax,nkptnr_loc))
allocate(docc(nmemax,nkptnr_loc))
call getmeidx(.false.,occsvnr)
call timer_stop(1)
if (wproc) then
  write(150,*)
  write(150,'("Maximum number of interband transitions: ",I5)')nmemax
  write(150,'("Done in ",F8.2," seconds")')timer(1,2)
endif
if (spinpol) then
  deallocate(spinor_ud)
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

call qname(ivq0m,qnm)
fname=trim(qnm)//"_me.hdf5"

if (root_cart((/1,1,0/))) then
  call h5fcreate_f(trim(fname),H5F_ACC_TRUNC_F,h5_root_id,ierr)
  call h5gcreate_f(h5_root_id,'parameters',h5_tmp_id,ierr)
  call h5gclose_f(h5_tmp_id,ierr)
  call h5gcreate_f(h5_root_id,'kpoints',h5_kpoints_id,ierr)
  do ik=1,nkptnr
    write(c8,'(I8.8)')ik
    call h5gcreate_f(h5_kpoints_id,c8,h5_kpoint_id,ierr)
    call h5gclose_f(h5_kpoint_id,ierr)
  enddo
  call h5gclose_f(h5_kpoints_id,ierr)
  call h5fclose_f(h5_root_id,ierr)
  call write_integer(nkptnr,1,trim(fname),'/parameters','nkptnr')
  call write_integer(nmemax,1,trim(fname),'/parameters','nmemax')
  call write_integer(lr_igq0,1,trim(fname),'/parameters','lr_igq0')
  call write_integer(gshme1,1,trim(fname),'/parameters','gshme1')
  call write_integer(gshme2,1,trim(fname),'/parameters','gshme2')
  call write_integer(gvecme1,1,trim(fname),'/parameters','gvecme1')
  call write_integer(gvecme2,1,trim(fname),'/parameters','gvecme2')
  call write_integer(ngvecme,1,trim(fname),'/parameters','ngvecme')
  call write_integer(nspinor,1,trim(fname),'/parameters','nspinor')
  call write_integer(spin_me,1,trim(fname),'/parameters','spin_me')
  call write_real8(vq0l,3,trim(fname),'/parameters','vq0l')
  call write_real8(vq0rl,3,trim(fname),'/parameters','vq0rl')
  call write_real8(vq0c,3,trim(fname),'/parameters','vq0c')
  call write_real8(vq0rc,3,trim(fname),'/parameters','vq0rc')
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

allocate(me(ngvecme,nmemax))
allocate(wfsvmt2(lmmaxvr,nrfmax,natmtot,nstsv,nspinor))
allocate(wfsvit2(ngkmax,nstsv,nspinor))
allocate(igkignr2(ngkmax))
allocate(evecfv2(nmatmax,nstfv,nspnfv))
allocate(evecsv2(nstsv,nstsv))

if (wproc) then
  write(150,*)
  write(150,'("Starting k-point loop")')
endif
do ikstep=1,nkptnrloc(0)
  call timer_reset(1)
  call timer_reset(2)
  call timer_reset(3)
  call timer_reset(4)
  call timer_reset(5) 
  if (wproc) then
    write(150,'("k-step ",I4," out of ",I4)')ikstep,nkptnrloc(0)
    call flushifc(150)
  endif
! transmit wave-functions
  call timer_start(1)
  call wfkq(ikstep,wfsvmtloc,wfsvitloc,ngknr,igkignr,wfsvmt2, &
    wfsvit2,ngknr2,igkignr2,evecfv2,evecsv2)
  call barrier(comm_cart)
  if (wproc) then
    write(150,'("  wave-functions are distributed")')
    call flushifc(150)
  endif
  call timer_stop(1)
! compute matrix elements  
  call timer_start(2)
  if (ikstep.le.nkptnrloc(mpi_x(1))) then
    call idxglob(nkptnr,mpi_dims(1),mpi_x(1)+1,ikstep,ik)
    if (.not.meoff) then
      me=dcmplx(0.d0,0.d0)
! calculate muffin-tin contribution for all combinations of n,n'    
      call timer_start(4)
      call zrhoftmt(nme(ikstep),ime(1,1,ikstep),               &
        wfsvmtloc(1,1,1,1,1,ikstep),wfsvmt2,ngumax,ngu,gu,igu, &
        me)
      call timer_stop(4)
! calculate interstitial contribution for all combinations of n,n'
      call timer_start(5)
      call zrhoftit(nme(ikstep),ime(1,1,ikstep),ngknr(ikstep), &
        ngknr2,igkignr(1,ikstep),igkignr2,idxkq(2,ik),         &
        wfsvitloc(1,1,1,ikstep),wfsvit2,me,evecfvloc(1,1,1,ikstep), &
        evecsvloc(1,1,ikstep),evecfv2,evecsv2,igfft1(1,ik))
      call timer_stop(5)
   else
      me=dcmplx(1.d0,0.d0)
    endif
  endif ! ikstep.le.nkptnrloc(mpi_x(1))
  call timer_stop(2)
! write matrix elements
  call timer_start(3)
  if (root_cart((/0,1,0/))) then
    do i=0,mpi_dims(1)-1
      do j=0,mpi_dims(3)-1
        if (i.eq.mpi_x(1).and.j.eq.mpi_x(3).and.ikstep.le.nkptnr_loc) then
          write(path,'("/kpoints/",I8.8)')ik
          call write_integer(idxkq(1,ik),1,trim(fname),trim(path),'kq')
          call write_integer(nme(ikstep),1,trim(fname),trim(path),'nme')
          call write_integer_array(ime(1,1,ikstep),2,(/3,nme(ikstep)/), &
            trim(fname),trim(path),'ime')
          call write_real8(docc(1,ikstep),nme(ikstep), &
            trim(fname),trim(path),'docc')
          call write_real8_array(me,3,(/2,ngvecme,nme(ikstep)/), &
            trim(fname),trim(path),'me')
        endif 
        call barrier(comm_cart_101)
      enddo
    enddo
  endif
  call timer_stop(3)
  if (wproc) then
    write(150,'("Time (seconds)")')
    write(150,'("  send/recv      : ",F8.2)')timer(1,2)
    write(150,'("  matrix elements")')
    write(150,'("    muffin-tins  : ",F8.2)')timer(4,2)
    write(150,'("    interstitial : ",F8.2)')timer(5,2)
    write(150,'("      total      : ",F8.2)')timer(2,2)
    write(150,'("    writing      : ",F8.2)')timer(3,2)
    write(150,'("  speed (me/sec) : ",F10.2)')ngvecme*nme(ikstep)/timer(2,2)
    call flushifc(150)
  endif
enddo !ikstep

deallocate(me)
deallocate(wfsvmt2)
deallocate(wfsvit2)
deallocate(igkignr2)
deallocate(evecfv2)
deallocate(evecsv2)

deallocate(nme)
deallocate(ime)
deallocate(docc)

deallocate(idxkq)
deallocate(vgq0c)
deallocate(gq0)
deallocate(tpgq0)
deallocate(sfacgq0)
deallocate(ylmgq0) 

deallocate(ngu)
deallocate(gu)
deallocate(igu)


if (wproc) then
  write(150,*)
  write(150,'("Done.")')
  call flushifc(150)
endif

return
end
#endif
