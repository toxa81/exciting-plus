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

! allocatable arrays
integer, allocatable :: igkignr2(:)
complex(8), allocatable :: wfsvmt2(:,:,:,:,:)
complex(8), allocatable :: wfsvit2(:,:,:)
complex(8), allocatable :: me(:,:)

integer i,j,ik,jk,ig,ikstep,ierr
integer ngknr2
real(8) vkq0l(3)
integer ivg1(3)
real(8) cpu0,cpu1,timeistl,timemt
real(8), allocatable :: uuj(:,:,:,:,:,:,:)
complex(4), allocatable :: gu(:,:,:)
integer, allocatable :: igu(:,:,:,:)
integer, allocatable :: ngu(:,:)
integer ngumax

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
enddo

! setup n,n' stuff
call getmeidx(.true.,occsvnr)
#ifdef _MPI_
if (root_cart((/0,1,1/))) then
  call mpi_allreduce(nmemax,i,1,MPI_INTEGER,MPI_MAX,comm_cart_100,ierr)
  nmemax=i
endif
call i_bcast_cart(comm_cart_011,nmemax,1)
#endif
allocate(nme(nkptnr_loc))
allocate(ime(3,nmemax,nkptnr_loc))
allocate(docc(nmemax,nkptnr_loc))
call getmeidx(.false.,occsvnr)
if (wproc) then
  write(150,*)
  write(150,'("Maximum number of interband transitions: ",I5)')nmemax
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

if (wproc) then
  write(150,*)
  write(150,'("Calculating radial integrals")')
  write(150,'("  maximum number of radial functions : ",I4)')nrfmax
endif
allocate(uuj(0:lmaxvr,0:lmaxvr,0:lmaxexp,nrfmax,nrfmax,natmtot,ngvecme))
call calc_uuj(uuj,lmaxexp,gq0)
if (wproc) then
  write(150,'("Done.")')
  call flushifc(150)
endif

call getgu(.true.,lmaxexp,uuj,ylmgq0,sfacgq0,ngumax,ngu,gu,igu)
allocate(ngu(natmtot,ngvecme))
allocate(gu(ngumax,natmtot,ngvecme))
allocate(igu(4,ngumax,natmtot,ngvecme))
call getgu(.false.,lmaxexp,uuj,ylmgq0,sfacgq0,ngumax,ngu,gu,igu)
if (wproc) then
  write(150,*)
  write(150,'("Maximum number of Gaunt-like coefficients : ",I8)')ngumax
  call flushifc(150)
endif
deallocate(uuj)

allocate(me(ngvecme,nmemax))
allocate(wfsvmt2(lmmaxvr,nrfmax,natmtot,nstsv,nspinor))
allocate(wfsvit2(ngkmax,nstsv,nspinor))
allocate(igkignr2(ngkmax))

if (wproc) then
  write(150,*)
  write(150,'("Starting k-point loop")')
endif
do ikstep=1,nkptnrloc(0)
  if (wproc) then
    write(150,'("k-step ",I4," out of ",I4)')ikstep,nkptnrloc(0)
    call flushifc(150)
  endif
! transmit wave-functions
  call wfkq(ikstep,wfsvmtloc,wfsvitloc,ngknr,igkignr,wfsvmt2, &
    wfsvit2,ngknr2,igkignr2)
  call barrier(comm_world)
  if (wproc) then
    write(150,'("  OK send and recieve")')
    call flushifc(150)
  endif

  if (ikstep.le.nkptnrloc(mpi_x(1))) then
    call idxglob(nkptnr,mpi_dims(1),mpi_x(1)+1,ikstep,ik)
    if (.not.meoff) then
      me=dcmplx(0.d0,0.d0)
! calculate interstitial contribution for all combinations of n,n'
      call cpu_time(cpu0)
      call zrhoftit(nme(ikstep),ime(1,1,ikstep),ngknr(ikstep), &
        ngknr2,igkignr(1,ikstep),igkignr2,idxkq(2,ik),         &
        wfsvitloc(1,1,1,ikstep),wfsvit2,me)
      call cpu_time(cpu1)
      timeistl=cpu1-cpu0
! calculate muffin-tin contribution for all combinations of n,n'    
      call cpu_time(cpu0)
      call zrhoftmt(nme(ikstep),ime(1,1,ikstep),               &
        wfsvmtloc(1,1,1,1,1,ikstep),wfsvmt2,ngumax,ngu,gu,igu, &
        me)
      call cpu_time(cpu1)
      timemt=cpu1-cpu0
      if (wproc) then
        write(150,'("  interstitial time (seconds) : ",F12.2)')timeistl
        write(150,'("    muffin-tin time (seconds) : ",F12.2)')timemt
        call flushifc(150)
      endif
    else
      me=dcmplx(1.d0,0.d0)
    endif
  endif ! (ikstep.le.nkptnrloc(iproc))
! write matrix elements
! todo: add timer count for writing
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
enddo !ikstep

deallocate(me)
deallocate(wfsvmt2)
deallocate(wfsvit2)
deallocate(igkignr2)

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


subroutine calc_uuj(uuj,lmaxexp,gq0)
use modmain
implicit none
! arguments
integer, intent(in) :: lmaxexp
real(8), intent(in) :: gq0(ngvecme)
real(8), intent(out) :: uuj(0:lmaxvr,0:lmaxvr,0:lmaxexp,nrfmax,nrfmax,natmtot,ngvecme)
! local variables
integer ia,is,ias,ig,l1,l2,l3,io1,io2,ir
real(8), allocatable :: jlgq0r(:,:,:,:)
real(8) fr(nrmtmax),gr(nrmtmax),cf(3,nrmtmax)
real(8) t1,jl(0:lmaxexp)

allocate(jlgq0r(nrmtmax,0:lmaxexp,nspecies,ngvecme))

! generate Bessel functions j_l(|G+q'|x)
do ig=1,ngvecme
  do is=1,nspecies
    do ir=1,nrmt(is)
      t1=gq0(ig)*spr(ir,is)
      call sbessel(lmaxexp,t1,jl)
      jlgq0r(ir,:,is,ig)=jl(:)
    enddo
  enddo
enddo

uuj=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ig=1,ngvecme
      do l1=0,lmaxvr
        do l2=0,lmaxvr
          do l3=0,lmaxexp
            do io1=1,nrfmax
              do io2=1,nrfmax
                do ir=1,nrmt(is)
                  fr(ir)=urf(ir,l1,io1,ias)*urf(ir,l2,io2,ias)*jlgq0r(ir,l3,is,ig)*(spr(ir,is)**2)
                enddo
                call fderiv(-1,nrmt(is),spr(1,is),fr,gr,cf)
                uuj(l1,l2,l3,io1,io2,ias,ig)=gr(nrmt(is))
              enddo !io2
            enddo !io1
          enddo !l3
        enddo !l2
      enddo !l1   
    enddo !ig
  enddo !ia
enddo !is

deallocate(jlgq0r)
return
end   

subroutine zrhoftit(nme0,ime0,ngknr1,ngknr2,igkignr1,igkignr2, &
  igkq,wfsvit1,wfsvit2,zrhofc0)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
integer, intent(in) :: nme0
integer, intent(in) :: ime0(3,nmemax)
integer, intent(in) :: ngknr1
integer, intent(in) :: ngknr2
integer, intent(in) :: igkq
integer, intent(in) :: igkignr1(ngkmax)
integer, intent(in) :: igkignr2(ngkmax)
complex(8), intent(in) :: wfsvit1(ngkmax,nstsv,nspinor)
complex(8), intent(in) :: wfsvit2(ngkmax,nstsv,nspinor)
complex(8), intent(inout) :: zrhofc0(ngvecme,nmemax)

complex(8), allocatable :: mit(:,:)
complex(8), allocatable :: a(:,:,:) 
complex(8), allocatable :: zrhofc_tmp(:,:)


integer is,ia,ias,ig,ig1,ig2,ist1,ist2,i,ispn,ispn2
integer idx0,bs,idx_g1,idx_g2
integer iv3g(3)
real(8) v1(3),v2(3),tp3g(2),len3g
complex(8) sfac3g(natmtot)

allocate(a(ngknr2,nstsv,nspinor))
allocate(zrhofc_tmp(ngvecme,nmemax))
zrhofc_tmp=dcmplx(0.d0,0.d0)

call idxbos(ngvecme,mpi_dims(2),mpi_x(2)+1,idx0,bs)
idx_g1=idx0+1
idx_g2=idx0+bs

allocate(mit(ngknr1,ngknr2))

do ig=idx_g1,idx_g2
  mit=dcmplx(0.d0,0.d0)
  do ig1=1,ngknr1
    do ig2=1,ngknr2
      ! G1-G2+G+K
      iv3g(:)=ivg(:,igkignr1(ig1))-ivg(:,igkignr2(ig2))+ivg(:,ig+gvecme1-1)+ivg(:,igkq)
      if (sum(abs(iv3g)).eq.0) mit(ig1,ig2)=dcmplx(1.d0,0.d0)
      v2(:)=1.d0*iv3g(:)
      call r3mv(bvec,v2,v1)
      call sphcrd(v1,len3g,tp3g)
      call gensfacgp(1,v1,1,sfac3g)
      do is=1,nspecies
        do ia=1,natoms(is)
	  ias=idxas(ia,is)
	  if (len3g.lt.1d-8) then
	    mit(ig1,ig2)=mit(ig1,ig2)-(fourpi/omega)*dconjg(sfac3g(ias))*(rmt(is)**3)/3.d0
	  else
	    mit(ig1,ig2)=mit(ig1,ig2)-(fourpi/omega)*dconjg(sfac3g(ias)) * &
	      (-(rmt(is)/len3g**2)*cos(len3g*rmt(is))+(1/len3g**3)*sin(len3g*rmt(is)))
	  endif
	enddo !ia
      enddo !is
    enddo
  enddo
  
  a=dcmplx(0.d0,0.d0)
  do ispn=1,nspinor
    do i=1,nstsv
      do ig2=1,ngknr2
        do ig1=1,ngknr1
          a(ig2,i,ispn)=a(ig2,i,ispn) + &
	    dconjg(wfsvit1(ig1,i,ispn))*mit(ig1,ig2)
        enddo
      enddo
    enddo
  enddo
  
  do ispn=1,nspinor
    if (lrtype.eq.0) then
      ispn2=ispn
    endif
    if (lrtype.eq.1) then
      ispn2=3-ispn
    endif
    do i=1,nme0
      ist1=ime0(1,i)
      ist2=ime0(2,i)
      do ig2=1,ngknr2
        zrhofc_tmp(ig,i)=zrhofc_tmp(ig,i)+wfsvit2(ig2,ist2,ispn2)*a(ig2,ist1,ispn)
      enddo
    enddo
  enddo
enddo !ig

if (mpi_dims(2).gt.1) then
  call d_reduce_cart(comm_cart_010,.false.,zrhofc_tmp,2*ngvecme*nmemax)
endif
zrhofc0=zrhofc0+zrhofc_tmp

deallocate(mit,a,zrhofc_tmp)
return
end
        
subroutine zrhoftmt(nme0,ime0,wfsvmt1,wfsvmt2,ngumax,ngu,gu,igu, &
  zrhofc0)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
integer, intent(in) :: nme0
integer, intent(in) :: ime0(3,nmemax)
integer, intent(in) :: ngumax
integer, intent(in) :: ngu(natmtot,ngvecme)
integer, intent(in) :: igu(4,ngumax,natmtot,ngvecme)
complex(4), intent(in) :: gu(ngumax,natmtot,ngvecme)
complex(8), intent(in) :: wfsvmt1(lmmaxvr,nrfmax,natmtot,nstsv,nspinor)
complex(8), intent(in) :: wfsvmt2(lmmaxvr,nrfmax,natmtot,nstsv,nspinor)
complex(8), intent(inout) :: zrhofc0(ngvecme,nmemax)
! local variables
integer ig,i,j,ist1,ist2,ias,io1,io2,lm1,lm2,ispn,ispn2
integer idx_g1,idx_g2,idx0,bs
complex(8) a1(lmmaxvr,nrfmax),a2(lmmaxvr,nrfmax)
complex(8), allocatable :: zrhofc_tmp(:,:)


allocate(zrhofc_tmp(ngvecme,nmemax))
zrhofc_tmp=dcmplx(0.d0,0.d0)

call idxbos(ngvecme,mpi_dims(2),mpi_x(2)+1,idx0,bs)
idx_g1=idx0+1
idx_g2=idx0+bs

do ig=idx_g1,idx_g2
  do ispn=1,nspinor
    if (lrtype.eq.0) then
      ispn2=ispn
    endif
    if (lrtype.eq.1) then
      ispn2=3-ispn
    endif
    do i=1,nme0
      ist1=ime0(1,i)
      ist2=ime0(2,i)
      do ias=1,natmtot
        a1=dconjg(wfsvmt1(:,:,ias,ist1,ispn))
        a2=wfsvmt2(:,:,ias,ist2,ispn2)
        do j=1,ngu(ias,ig)
          lm1=igu(1,j,ias,ig)
          lm2=igu(2,j,ias,ig)
          io1=igu(3,j,ias,ig)
          io2=igu(4,j,ias,ig)
          zrhofc_tmp(ig,i)=zrhofc_tmp(ig,i)+a1(lm1,io1)*a2(lm2,io2)*gu(j,ias,ig)
        enddo
      enddo !ias 
    enddo !i
  enddo
enddo !ig    

if (mpi_dims(2).gt.1) then
  call d_reduce_cart(comm_cart_010,.false.,zrhofc_tmp,2*ngvecme*nmemax)
endif
zrhofc0=zrhofc0+zrhofc_tmp

deallocate(zrhofc_tmp)

return
end


subroutine getmeidx(req,occsvnr)
use modmain
implicit none
! arguments
logical, intent(in) :: req
real(8), intent(in) :: occsvnr(nstsv,nkptnr)

integer band1,band2
integer i,ik,jk,istfv1,istfv2,ispn1,ispn2,ist1,ist2,ikloc
logical laddme,ldocc
real(8) d1
integer, external :: iknrglob2

if (bndme1.eq.-1) then
  band1=1
  band2=nstfv
else
  band1=bndme1
  band2=bndme2
endif
if (req.and.wproc) then
  write(150,*)
  write(150,'("Band interval: ",2I4)')band1,band2
endif
if (req) nmemax=0
do ikloc=1,nkptnr_loc
  ik=iknrglob2(ikloc,mpi_x(1))
  jk=idxkq(1,ik)
  i=0
  do ispn1=1,nspinor
    do ispn2=1,nspinor
      do istfv1=band1,band2
      do istfv2=band1,band2
        ist1=istfv1+(ispn1-1)*nstfv
        ist2=istfv2+(ispn2-1)*nstfv
        d1=occsvnr(ist1,ik)-occsvnr(ist2,jk)
        ldocc=abs(d1).gt.1d-10
        laddme=.false.
! possible candidate for charge response
        if (ispn1.eq.ispn2.and.lrtype.eq.0) then
          if ((ispn1.eq.spin_me.or.spin_me.eq.3).and.ldocc) laddme=.true.
        endif
! for magnetic response
        if (ispn1.ne.ispn2.and.lrtype.eq.1) then
          if ((ispn1.eq.spin_me.or.spin_me.eq.3).and.ldocc) laddme=.true.
        endif
        if (laddme) then
          i=i+1
          if (.not.req) then
            ime(1,i,ikloc)=ist1
            ime(2,i,ikloc)=ist2
            ime(3,i,ikloc)=ispn1
            docc(i,ikloc)=d1
          endif
        endif
      enddo !istfv2
      enddo !istfv1
    enddo !ispn2
  enddo !ispn1
  if (.not.req) nme(ikloc)=i
  if (req) nmemax=max(nmemax,i)
enddo !ikloc

return
end

subroutine getgu(req,lmaxexp,uuj,ylmgq0,sfacgq0,ngumax,ngu,gu,igu)
use modmain
implicit none
! arguments
logical, intent(in) :: req
integer, intent(in) :: lmaxexp
real(8), intent(in) :: uuj(0:lmaxvr,0:lmaxvr,0:lmaxexp,nrfmax,nrfmax,natmtot,ngvecme)
complex(8), intent(in) :: ylmgq0((lmaxexp+1)**2,ngvecme)
complex(8), intent(in) :: sfacgq0(ngvecme,natmtot)
integer, intent(inout) :: ngumax
integer, intent(out) :: ngu(natmtot,ngvecme)
complex(4), intent(out) :: gu(ngumax,natmtot,ngvecme)
integer, intent(out) :: igu(4,ngumax,natmtot,ngvecme)

integer ig,ias,i,io1,io2,l1,m1,lm1,l2,m2,lm2,l3,m3,lm3
real(8) t1
real(8), external :: gaunt

if (req) ngumax=0 

do ig=1,ngvecme
  do ias=1,natmtot
    i=0
    do io1=1,nrfmax
      do io2=1,nrfmax
        do l1=0,lmaxvr
        do m1=-l1,l1 
          lm1=idxlm(l1,m1)
          do l2=0,lmaxvr
          do m2=-l2,l2
            lm2=idxlm(l2,m2)
            do l3=0,lmaxexp
            do m3=-l3,l3
              lm3=idxlm(l3,m3)
	      t1=gaunt(l2,l1,l3,m2,m1,m3)*uuj(l1,l2,l3,io1,io2,ias,ig)
              if (abs(t1).gt.1d-10) then
                i=i+1
	        if (.not.req) then
                  gu(i,ias,ig)=t1*ylmgq0(lm3,ig)*dconjg(zi**l3)*fourpi*dconjg(sfacgq0(ig,ias))
                  igu(1,i,ias,ig)=lm1
                  igu(2,i,ias,ig)=lm2
                  igu(3,i,ias,ig)=io1
                  igu(4,i,ias,ig)=io2
	        endif
              endif
            enddo
            enddo
          enddo
          enddo
        enddo
        enddo
      enddo
    enddo
    if (.not.req) ngu(ias,ig)=i
    if (req) ngumax=max(ngumax,i)
  enddo
enddo    

return
end

subroutine writegshells
use modmain
implicit none

integer ngsh
integer ,allocatable :: igishell(:)
integer ,allocatable :: ishellng(:,:)

integer ig,ish

allocate(igishell(ngvec))
allocate(ishellng(ngvec,2))

call getgshells(ngsh,igishell,ishellng)
open(170,file='GSHELLS.OUT',form='formatted',status='replace')
write(170,*)
write(170,'("  G-vec.       lat.coord.      length(a.u.)   shell ")')
write(170,'(" ---------------------------------------------------")')
do ig=1,ishellng(ngsh,2)
  write(170,'(2X,I4,4X,3I5,4X,F12.6,5x,I4)')ig,ivg(:,ig),gc(ig),igishell(ig)
enddo
write(170,*)
write(170,'("  G-shell    num. G-vec in shell      total num. G-vec. ")')
write(170,'(" -------------------------------------------------------")')
do ish=1,ngsh
  write(170,'(2X,I4,14X,I4,18X,I4)')ish,ishellng(ish,1),ishellng(ish,2)
enddo
close(170)

deallocate(igishell,ishellng)

return
end


subroutine write_real8(a,n,fname,path,nm)
use hdf5
implicit none
integer, intent(in) :: n
real(8), intent(in) :: a(n)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm


integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr
integer(HSIZE_T), dimension(1) :: dims

dims(1)=n
call h5fopen_f(fname,H5F_ACC_RDWR_F,h5_root_id,ierr)
call h5screate_simple_f(1,dims,dataspace_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dcreate_f(group_id,nm,H5T_NATIVE_DOUBLE,dataspace_id, &
  dataset_id,ierr)
call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,a,dims,ierr)
call h5gclose_f(group_id,ierr)
call h5sclose_f(dataspace_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end

subroutine read_real8(a,n,fname,path,nm)
use hdf5
implicit none
integer, intent(in) :: n
real(8), intent(out) :: a(n)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm

integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr
integer(HSIZE_T), dimension(1) :: dims

dims(1)=n
call h5fopen_f(fname,H5F_ACC_RDONLY_F,h5_root_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dopen_f(group_id,nm,dataset_id,ierr)
call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,a,dims,ierr)
call h5gclose_f(group_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end


subroutine write_real8_array(a,ndims,dims,fname,path,nm)
use hdf5
implicit none
integer, intent(in) :: ndims
integer, intent(in) :: dims(ndims)
real(8), intent(in) :: a(*)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm


integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr,i
integer(HSIZE_T), dimension(ndims) :: h_dims

do i=1,ndims
  h_dims(i)=dims(i)
enddo
call h5fopen_f(fname,H5F_ACC_RDWR_F,h5_root_id,ierr)
call h5screate_simple_f(ndims,h_dims,dataspace_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dcreate_f(group_id,nm,H5T_NATIVE_DOUBLE,dataspace_id, &
  dataset_id,ierr)
call h5dwrite_f(dataset_id,H5T_NATIVE_DOUBLE,a,h_dims,ierr)
call h5gclose_f(group_id,ierr)
call h5sclose_f(dataspace_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end

subroutine read_real8_array(a,ndims,dims,fname,path,nm)
use hdf5
implicit none
integer, intent(in) :: ndims
integer, intent(in) :: dims(ndims)
real(8), intent(out) :: a(*)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm


integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr,i
integer(HSIZE_T), dimension(ndims) :: h_dims

do i=1,ndims
  h_dims(i)=dims(i)
enddo

call h5fopen_f(fname,H5F_ACC_RDONLY_F,h5_root_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dopen_f(group_id,nm,dataset_id,ierr)
call h5dread_f(dataset_id,H5T_NATIVE_DOUBLE,a,h_dims,ierr)
call h5gclose_f(group_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end


subroutine write_integer_array(a,ndims,dims,fname,path,nm)
use hdf5
implicit none
integer, intent(in) :: ndims
integer, intent(in) :: dims(ndims)
integer, intent(in) :: a(*)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm


integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr,i
integer(HSIZE_T), dimension(ndims) :: h_dims

do i=1,ndims
  h_dims(i)=dims(i)
enddo
call h5fopen_f(fname,H5F_ACC_RDWR_F,h5_root_id,ierr)
call h5screate_simple_f(ndims,h_dims,dataspace_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dcreate_f(group_id,nm,H5T_NATIVE_INTEGER,dataspace_id, &
  dataset_id,ierr)
call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,a,h_dims,ierr)
call h5gclose_f(group_id,ierr)
call h5sclose_f(dataspace_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end


subroutine read_integer_array(a,ndims,dims,fname,path,nm)
use hdf5
implicit none
integer, intent(in) :: ndims
integer, intent(in) :: dims(ndims)
integer, intent(out) :: a(*)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm


integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr,i
integer(HSIZE_T), dimension(ndims) :: h_dims

do i=1,ndims
  h_dims(i)=dims(i)
enddo
call h5fopen_f(fname,H5F_ACC_RDONLY_F,h5_root_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dopen_f(group_id,nm,dataset_id,ierr)
call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,a,h_dims,ierr)
call h5gclose_f(group_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end

subroutine write_integer(a,n,fname,path,nm)
use hdf5
implicit none
integer, intent(in) :: n
integer, intent(in) :: a(n)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm


integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr
integer(HSIZE_T), dimension(1) :: dims

dims(1)=n
call h5fopen_f(fname,H5F_ACC_RDWR_F,h5_root_id,ierr)
call h5screate_simple_f(1,dims,dataspace_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dcreate_f(group_id,nm,H5T_NATIVE_INTEGER,dataspace_id, &
  dataset_id,ierr)
call h5dwrite_f(dataset_id,H5T_NATIVE_INTEGER,a,dims,ierr)
call h5gclose_f(group_id,ierr)
call h5sclose_f(dataspace_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end

subroutine read_integer(a,n,fname,path,nm)
use hdf5
implicit none
integer, intent(in) :: n
integer, intent(out) :: a(n)
character(*), intent(in) :: fname
character(*), intent(in) :: path
character(*), intent(in) :: nm

integer(hid_t) h5_root_id,dataspace_id,dataset_id,group_id
integer ierr
integer(HSIZE_T), dimension(1) :: dims

dims(1)=n
call h5fopen_f(fname,H5F_ACC_RDONLY_F,h5_root_id,ierr)
call h5gopen_f(h5_root_id,path,group_id,ierr)
call h5dopen_f(group_id,nm,dataset_id,ierr)
call h5dread_f(dataset_id,H5T_NATIVE_INTEGER,a,dims,ierr)
call h5gclose_f(group_id,ierr)
call h5dclose_f(dataset_id,ierr)
call h5fclose_f(h5_root_id,ierr)
return
end














subroutine wfkq(ikstep,wfsvmtloc,wfsvitloc,ngknr,igkignr,wfsvmt2, &
  wfsvit2,ngknr2,igkignr2)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
integer, intent(in) :: ikstep
complex(8), intent(in) :: wfsvmtloc(lmmaxvr,nrfmax,natmtot,nstsv,nspinor,*)
complex(8), intent(in) :: wfsvitloc(ngkmax,nstsv,nspinor,*)
integer, intent(in) :: ngknr(*)
integer, intent(in) :: igkignr(ngkmax,*)
complex(8), intent(out) :: wfsvmt2(lmmaxvr,nrfmax,natmtot,nstsv,nspinor)
complex(8), intent(out) :: wfsvit2(ngkmax,nstsv,nspinor)
integer, intent(out) :: ngknr2
integer, intent(out) :: igkignr2(ngkmax)


integer idx0,bs,ip
logical lsend,lrecv
integer i,ik,jk,ikloc
integer rank,tag,ierr,req
integer, allocatable :: stat(:)

#ifdef _MPI_
! each proc in k-group knows that it needs wave-function at jk=idxkq(1,ik) point
!
! the distribution of k-points could look like this
!                p0          p1          p2
!          +-----------+-----------+-----------+
! ikstep=1 | ik=1 jk=3 | ik=4 jk=2 | ik=7 jk=5 |
! ikstep=2 | ik=2 jk=4 | ik=5 jk=7 | ik=8 jk=6 |
! ikstep=3 | ik=3 jk=1 | ik=6 jk=8 |  -        |
!          +-----------+-----------+-----------+
allocate(stat(MPI_STATUS_SIZE))
do i=1,mpi_dims(1)
  call idxbos(nkptnr,mpi_dims(1),i,idx0,bs)
  if (ikstep.le.bs) then
! for step ikstep process i handles k-point ik
    call idxglob(nkptnr,mpi_dims(1),i,ikstep,ik)
! for step ikstep process i requires k-point jk
    jk=idxkq(1,ik)
! find the process which stores the k-point jk and it's local index
    call idxloc(nkptnr,mpi_dims(1),jk,ip,ikloc)
! send/recv when ip.ne.i (when k-points ik and jk are on different procs)
    if (mpi_x(1).eq.ip-1.and.ip.ne.i) then
! destination proc is i-1
      call mpi_cart_rank(comm_cart_100,(/i-1/),rank,ierr) 
! send muffin-tin part
      tag=(ikstep*nproc+i-1)*10
      call mpi_isend(wfsvmtloc(1,1,1,1,1,ikloc),                     &
        lmmaxvr*nrfmax*natmtot*nstsv*nspinor,                        & 
	    MPI_DOUBLE_COMPLEX,rank,tag,comm_cart_100,req,ierr)
! send interstitial part
      tag=tag+1
      call mpi_isend(wfsvitloc(1,1,1,ikloc),ngkmax*nstsv*nspinor,    & 
        MPI_DOUBLE_COMPLEX,rank,tag,comm_cart_100,req,ierr)
! send ngknr      
      tag=tag+1
      call mpi_isend(ngknr(ikloc),1,MPI_INTEGER,rank,tag,comm_cart_100, &
        req,ierr)
! send igkignr
      tag=tag+1
      call mpi_isend(igkignr(1,ikloc),ngkmax,MPI_INTEGER,rank,tag,   &
        comm_cart_100,req,ierr)  
    endif
    if (mpi_x(1).eq.i-1.and.ip.ne.i) then
! source proc is ip-1
      call mpi_cart_rank(comm_cart_100,(/ip-1/),rank,ierr)
      tag=(ikstep*nproc+mpi_x(1))*10
      call mpi_recv(wfsvmt2,lmmaxvr*nrfmax*natmtot*nstsv*nspinor,    &
        MPI_DOUBLE_COMPLEX,rank,tag,comm_cart_100,stat,ierr)
      tag=tag+1
      call mpi_recv(wfsvit2,ngkmax*nstsv*nspinor,MPI_DOUBLE_COMPLEX, &
        rank,tag,comm_cart_100,stat,ierr)
      tag=tag+1
      call mpi_recv(ngknr2,1,MPI_INTEGER,rank,tag,comm_cart_100,stat,   &
        ierr)
      tag=tag+1
      call mpi_recv(igkignr2,ngkmax,MPI_INTEGER,rank,tag,comm_cart_100, &
        stat,ierr)
    endif
    if (mpi_x(1).eq.i-1.and.ip.eq.i) then
      wfsvmt2(:,:,:,:,:)=wfsvmtloc(:,:,:,:,:,ikloc)
      wfsvit2(:,:,:)=wfsvitloc(:,:,:,ikloc)
      ngknr2=ngknr(ikloc)
      igkignr2(:)=igkignr(:,ikloc)
    endif
  endif
enddo
deallocate(stat)
#else
  jk=idxkq(1,ikstep)
  wfsvmt2(:,:,:,:,:)=wfsvmtloc(:,:,:,:,:,jk)
  wfsvit2(:,:,:)=wfsvitloc(:,:,:,jk)
  ngknr2=ngknr(jk)
  igkignr2(:)=igkignr(:,jk)
#endif
return
end




subroutine qname(ivq,name)
implicit none
integer, intent(in) :: ivq(3)
character*(*), intent(out) :: name
character*6 fmt1,str1
integer i
name='q_'
do i=1,3
  if (ivq(i).lt.0) then 
    fmt1='(I4.3)'
  else
    fmt1='(I3.3)'
  endif
  write(str1,fmt1)ivq(i)
  name=trim(name)//trim(adjustl(str1))
  if (i.lt.3) name=trim(name)//"_"
enddo
return
end
