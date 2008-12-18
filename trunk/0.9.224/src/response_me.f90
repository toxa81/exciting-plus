subroutine response_me(ivq0m,gvecme1,gvecme2,ngvecme,wfsvmtloc, &
  wfsvitloc,ngknr,igkignr,occsvnr)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
! q-vector in k-mesh coordinates
integer, intent(in) :: ivq0m(3)
integer, intent(in) :: gvecme1
integer, intent(in) :: gvecme2
! number of G-vectors for matrix elements calculation
integer, intent(in) :: ngvecme
complex(8), intent(in) :: wfsvmtloc(lmmaxvr,nrfmax,natmtot,nstsv,nspinor,*)
complex(8), intent(in) :: wfsvitloc(nmatmax,nstsv,nspinor,*)
integer, intent(in) :: ngknr(*)
integer, intent(in) :: igkignr(ngkmax,*)
real(8), intent(in) :: occsvnr(nstsv,nkptnr)
 
! q-vector in lattice coordinates
real(8) vq0l(3)
! q-vector in Cartesian coordinates
real(8) vq0c(3)
! reduced q-vector in lattice coordinates
real(8) vq0rl(3)
! reduced q-vector in Cartesian coordinates
real(8) vq0rc(3)
! G-vector which brings q to first BZ
integer vgq0l(3)
! index of G-vector which brings q to first BZ
integer igq0
! array for k and k+q stuff
!  1-st index: index of k-point in BZ
!  2-nd index: 1: index of k'=k+q-K
!              2: index of K-vector which brings k+q to first BZ
integer, allocatable :: ikq(:,:)
! number of n,n' combinations of band indexes for each k-point
integer, allocatable :: num_nnp(:)
! maximum num_nnp over all k-points 
integer max_num_nnp
! pair of n,n' band indexes for each k-point
integer, allocatable :: nnp(:,:,:)
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
real(8), allocatable :: docc(:,:)
complex(8), allocatable :: zrhofc(:,:,:)
complex(8), allocatable :: wfsvmt2(:,:,:,:,:)
complex(8), allocatable :: wfsvit2(:,:,:)

integer i,j,i1,ik,jk,ig,is,ikstep,ist1,ist2,ispn,ikloc,l,ia,ias,istfv
integer ngknr2
real(8) vkq0l(3),t1,jl(0:lmaxvr)
integer ivg1(3),ivg2(3)
real(8) cpu0,cpu1,timeistl,timemt
real(8), allocatable :: uuj(:,:,:,:,:,:,:)
complex(4), allocatable :: gu(:,:,:)
integer, allocatable :: igu(:,:,:,:)
integer, allocatable :: ngu(:,:)
integer ngumax

integer lmaxexp
integer lmmaxexp

! for parallel execution
integer, allocatable :: ikptiprocnr(:)
integer, allocatable :: isend(:,:,:)
integer tag,req,ierr
integer, allocatable :: stat(:)
character*100 :: fname

! external functions
real(8), external :: r3taxi
complex(8), external :: zfint

lmaxexp=lmaxvr
lmmaxexp=(lmaxexp+1)**2

if (iproc.eq.0) then
  write(150,*)
  write(150,'("Calculation of matrix elements <n,k|e^{-i(G+q)x}|n'',k+q>")')
  if (spinpol) then
    write(150,*)
    write(150,'("Spin-polarized calculation")')
    if (lrtype.eq.0) then
      if (spin_me.eq.1) write(150,'(" calculation of matrix elements for spin up")')
      if (spin_me.eq.2) write(150,'(" calculation of matrix elements for spin dn")')
      if (spin_me.eq.3) write(150,'(" calculation of matrix elements for both spins")')
    endif
  endif
endif

allocate(ikptiprocnr(nkptnr))
do i=0,nproc-1
  ikptiprocnr(ikptnrloc(i,1):ikptnrloc(i,2))=i
enddo

allocate(ikq(nkptnr,2))
allocate(vgq0c(3,ngvecme))
allocate(gq0(ngvecme))
allocate(tpgq0(2,ngvecme))
allocate(sfacgq0(ngvecme,natmtot))
allocate(ylmgq0(lmmaxexp,ngvecme)) 


! q-vector in lattice coordinates
do i=1,3
  vq0l(i)=1.d0*ivq0m(i)/ngridk(i)+1d-12
enddo

! find G-vector which brings q0 to first BZ
vgq0l(:)=floor(vq0l(:))

! reduce q0 vector to first BZ
vq0rl(:)=vq0l(:)-vgq0l(:)

! check if we have enough G-shells to bring q-vector back to first BZ
do ig=1,ngvecme
  if (sum(abs(vgq0l(:)-ivg(:,ig+gvecme1-1))).eq.0) then
    igq0=ig+gvecme1-1
    goto 20
  endif
enddo
write(*,*)
write(*,'("Bug(response_me): no G-vector to reduce q-vector to first BZ")')
write(*,*)
call pstop
20 continue

! get Cartesian coordinates of q-vector and reduced q-vector
call r3mv(bvec,vq0l,vq0c)
call r3mv(bvec,vq0rl,vq0rc)
  
if (iproc.eq.0) then
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
    & I4)')igq0
  write(150,'("G-vector (lat.coord.)                        : ",&
    & 3I4)')ivg(:,igq0)
  write(150,'("reduced q-vector (lat.coord.)                : ",&
    & 3G18.10)')vq0rl
  write(150,'("reduced q-vector (Cart.coord.) [a.u.]        : ",&
    & 3G18.10)')vq0rc
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
      ikq(ik,1)=jk
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
  ikq(ik,2)=ivgig(ivg1(1),ivg1(2),ivg1(3))
enddo

! setup n,n' stuff
call getnnp(.true.,ikq,occsvnr,max_num_nnp,num_nnp,nnp,docc)
#ifdef _MPI_
call mpi_allreduce(max_num_nnp,i,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
max_num_nnp=i
#endif  
allocate(num_nnp(nkptnrloc(iproc)))
allocate(nnp(max_num_nnp,3,nkptnrloc(iproc)))
allocate(docc(max_num_nnp,nkptnrloc(iproc)))
call getnnp(.false.,ikq,occsvnr,max_num_nnp,num_nnp,nnp,docc)
if (iproc.eq.0) then
  write(150,*)
  write(150,'("Maximum number of interband transitions: ",I5)')max_num_nnp
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

write(fname,'("ZRHOFC[",I4.3,",",I4.3,",",I4.3,"].OUT")') &
  ivq0m(1),ivq0m(2),ivq0m(3)

if (iproc.eq.0) then
  open(160,file=trim(fname),form='unformatted',status='replace')
  write(160)nkptnr,max_num_nnp,igq0
  write(160)gshme1,gshme2,gvecme1,gvecme2,ngvecme
  write(160)nspinor,spin_me
  write(160)vq0l(1:3)
  write(160)vq0rl(1:3)
  write(160)vq0c(1:3)
  write(160)vq0rc(1:3)
  do ik=1,nkptnr
    write(160)ikq(ik,1)
  enddo
  close(160)
endif

if (iproc.eq.0) then
  write(150,*)
  write(150,'("Calculating radial integrals")')
  write(150,'("  maximum number of radial functions : ",I4)')nrfmax
endif
allocate(uuj(0:lmaxvr,0:lmaxvr,0:lmaxexp,nrfmax,nrfmax,natmtot,ngvecme))
call calc_uuj(uuj,lmaxexp,gq0,ngvecme)
if (iproc.eq.0) then
  write(150,'("Done.")')
  call flushifc(150)
endif

call getgu(.true.,ngvecme,lmaxexp,uuj,ylmgq0,sfacgq0,ngumax,ngu,gu,igu)
allocate(ngu(natmtot,ngvecme))
allocate(gu(ngumax,natmtot,ngvecme))
allocate(igu(4,ngumax,natmtot,ngvecme))
call getgu(.false.,ngvecme,lmaxexp,uuj,ylmgq0,sfacgq0,ngumax,ngu,gu,igu)
if (iproc.eq.0) then
  write(150,*)
  write(150,'("Maximum number of Gaunt-like coefficients : ",I8)')ngumax
  call flushifc(150)
endif
deallocate(uuj)

! different implementation for parallel and serial execution
#ifdef _MPI_

allocate(wfsvmt2(lmmaxvr,nrfmax,natmtot,nstsv,nspinor))
allocate(wfsvit2(nmatmax,nstsv,nspinor))
allocate(igkignr2(ngkmax))
allocate(stat(MPI_STATUS_SIZE))
allocate(zrhofc(ngvecme,max_num_nnp,nkptnrloc(iproc)))


! do ikloc=1,nkptnrloc(iproc)
! hack to switch off l-channel
!  do is=1,nspecies
!    if (trim(adjustl(spsymb(is))).eq.'O') then
!      do ia=1,natoms(is)
!        ias=idxas(ia,is)
!	  acoeffloc(2:4,:,ias,:,29:55,ikloc)=dcmplx(0.d0,0.d0)
!      enddo
!    endif
!  enddo
!  ias=7
!  acoeffloc(2:4,:,ias,:,29:55,ikloc)=dcmplx(0.d0,0.d0)
!  ias=8
!  acoeffloc(2:4,:,ias,:,29:55,ikloc)=dcmplx(0.d0,0.d0)
! check p-d transition
!  ias=11
!  acoeffloc(2:4,:,ias,:,29:55,ikloc)=dcmplx(0.d0,0.d0)
!  ias=12
!  acoeffloc(2:4,:,ias,:,29:55,ikloc)=dcmplx(0.d0,0.d0)
!  ias=5
!  acoeffloc(5:9,:,ias,:,56,ikloc)=dcmplx(0.d0,0.d0)
! check d-p transition
!  ias=11
!  acoeffloc(2:4,:,ias,:,56,ikloc)=dcmplx(0.d0,0.d0)
!  ias=12
!  acoeffloc(2:4,:,ias,:,56,ikloc)=dcmplx(0.d0,0.d0)
!  ias=6
!  acoeffloc(5:9,:,ias,:,29:55,ikloc)=dcmplx(0.d0,0.d0)
! enddo

! find indexes of k-points to send and receive
allocate(isend(nkptnrloc(0),0:nproc-1,2))
isend=-1
do ikstep=1,nkptnrloc(0)
  do i=0,nproc-1
    ik=ikptnrloc(i,1)+ikstep-1
    if (ikstep.le.nkptnrloc(i)) then
! for non reduced k' point find the proc and local index
      isend(ikstep,i,1)=ikptiprocnr(ikq(ik,1))
      isend(ikstep,i,2)=ikq(ik,1)-ikptnrloc(isend(ikstep,i,1),1)+1
    endif
  enddo
enddo

if (iproc.eq.0) then
  write(150,*)
  write(150,'("Starting k-point loop")')
endif
zrhofc=dcmplx(0.d0,0.d0)
do ikstep=1,nkptnrloc(0)
  if (iproc.eq.0) then
    write(150,'("k-step ",I4," out of ",I4)')ikstep,nkptnrloc(0)
    call flushifc(150)
  endif
! find the i-th proc to which the current iproc should send data
  do i=0,nproc-1
    if (isend(ikstep,i,1).eq.iproc.and.iproc.ne.i) then
      tag=(ikstep*nproc+i)*10
      ik=isend(ikstep,i,2)
      call mpi_isend(wfsvitloc(1,1,1,ik),nmatmax*nstsv*nspinor,        &
        MPI_DOUBLE_COMPLEX,i,tag,MPI_COMM_WORLD,req,ierr)
      tag=tag+1
      call mpi_isend(wfsvmtloc(1,1,1,1,1,ik),                  &
        lmmaxvr*nrfmax*natmtot*nstsv*nspinor,                        &
	MPI_DOUBLE_COMPLEX,i,tag,MPI_COMM_WORLD,req,ierr)
      tag=tag+1
      call mpi_isend(ngknr(ik),1,MPI_INTEGER,i,tag,          &
        MPI_COMM_WORLD,req,ierr)
      tag=tag+1
      call mpi_isend(igkignr(1,ik),ngkmax,MPI_INTEGER,i,tag, &
        MPI_COMM_WORLD,req,ierr)  
    endif
  enddo
! receive data
  if (isend(ikstep,iproc,1).ne.-1) then
    if (isend(ikstep,iproc,1).ne.iproc) then
      tag=(ikstep*nproc+iproc)*10
      call mpi_recv(wfsvit2,nmatmax*nstsv*nspinor,MPI_DOUBLE_COMPLEX,  &
        isend(ikstep,iproc,1),tag,MPI_COMM_WORLD,stat,ierr)
      tag=tag+1
      call mpi_recv(wfsvmt2,lmmaxvr*nrfmax*natmtot*nstsv*nspinor,      &
        MPI_DOUBLE_COMPLEX,isend(ikstep,iproc,1),tag,          &
        MPI_COMM_WORLD,stat,ierr)
      tag=tag+1
      call mpi_recv(ngknr2,1,MPI_INTEGER,isend(ikstep,iproc,1), &
        tag,MPI_COMM_WORLD,stat,ierr)
      tag=tag+1
      call mpi_recv(igkignr2,ngkmax,MPI_INTEGER,isend(ikstep,iproc,1), &
        tag,MPI_COMM_WORLD,stat,ierr)
    else
      ik=isend(ikstep,iproc,2)
      wfsvit2(:,:,:)=wfsvitloc(:,:,:,ik)
      wfsvmt2(:,:,:,:,:)=wfsvmtloc(:,:,:,:,:,ik)
      ngknr2=ngknr(ik)
      igkignr2(:)=igkignr(:,ik)
    endif
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  if (iproc.eq.0) then
    write(150,'("  OK send and recieve")')
    call flushifc(150)
  endif
  
  if (ikstep.le.nkptnrloc(iproc)) then
    ik=ikptnrloc(iproc,1)+ikstep-1
    jk=ikq(ik,1)
    
    if (.not.meoff) then
! calculate interstitial contribution for all combinations of n,n'
      call cpu_time(cpu0)
      call zrhoftit(ngvecme,gvecme1,max_num_nnp,num_nnp(ikstep), &
        nnp(1,1,ikstep),ngknr(ikstep),ngknr2,igkignr(1,ikstep),  &
        igkignr2,ikq(ik,2),wfsvitloc(1,1,1,ikstep),wfsvit2,zrhofc(1,1,ikstep))
      call cpu_time(cpu1)
      timeistl=cpu1-cpu0

! calculate muffin-tin contribution for all combinations of n,n'    
      call cpu_time(cpu0)
      call zrhoftmt(ngvecme,max_num_nnp,num_nnp(ikstep),nnp(1,1,ikstep), &
        wfsvmtloc(1,1,1,1,1,ikstep),wfsvmt2,ngumax,ngu,gu,igu,zrhofc(1,1,ikstep))
      call cpu_time(cpu1)
      timemt=cpu1-cpu0
  
      if (iproc.eq.0) then
        write(150,'("  interstitial time (seconds) : ",F12.2)')timeistl
        write(150,'("    muffin-tin time (seconds) : ",F12.2)')timemt
        call flushifc(150)
      endif
    else
      zrhofc(:,:,ikstep)=dcmplx(1.d0,0.d0)
    endif
  endif ! (ikstep.le.nkptnrloc(iproc))
  
enddo !ikstep

call mpi_barrier(MPI_COMM_WORLD,ierr)

do i=0,nproc-1
  if (i.eq.iproc) then
    open(160,file=trim(fname),form='unformatted',status='old',position='append')
    do ikstep=1,nkptnrloc(iproc)
      ik=ikptnrloc(iproc,1)+ikstep-1
      write(160)ik,ikq(ik,1)
      write(160)num_nnp(ikstep)
      write(160)nnp(1:num_nnp(ikstep),1:3,ikstep)
      write(160)docc(1:num_nnp(ikstep),ikstep)
      write(160)zrhofc(1:ngvecme,1:num_nnp(ikstep),ikstep)
    enddo !ikstep
    close(160)
  endif !i.eq.iproc
  call mpi_barrier(MPI_COMM_WORLD,ierr)
enddo 
deallocate(wfsvmt2)
deallocate(wfsvit2)
deallocate(igkignr2)
deallocate(stat)
deallocate(zrhofc)
  
#else

allocate(zrhofc(ngvecme,max_num_nnp,1))

open(160,file=trim(fname),form='unformatted',status='old',position='append')

write(150,*)
do ik=1,nkptnr
  write(150,'("k-point ",I4," out of ",I4)')ik,nkptnr
  call flushifc(150)
  
  jk=ikq(ik,1)
  
  write(160)ik,jk
  write(160)num_nnp(ik)
  write(160)nnp(1:num_nnp(ik),1:3,ik)
  write(160)docc(1:num_nnp(ik),ik)
  
  zrhofc=dcmplx(0.d0,0.d0)
! calculate interstitial contribution for all combinations of n,n'
  call cpu_time(cpu0)
  call zrhoftit(ngvecme,gvecme1,max_num_nnp,num_nnp(ik),nnp(1,1,ik), &
    ngknr(ik),ngknr(jk),igkignr(1,ik),igkignr(1,jk),ikq(ik,2),       &
    wfsvitloc(1,1,1,ik),wfsvitloc(1,1,1,jk),zrhofc)
  call cpu_time(cpu1)
  timeistl=cpu1-cpu0

! calculate muffin-tin contribution for all combinations of n,n'    
  call cpu_time(cpu0)
  call zrhoftmt(ngvecme,max_num_nnp,num_nnp(ik),nnp(1,1,ik), &
    wfsvmtloc(1,1,1,1,1,ik),wfsvmtloc(1,1,1,1,1,jk),ngumax,ngu,gu,igu,zrhofc)
  call cpu_time(cpu1)
  timemt=cpu1-cpu0
  
  if (iproc.eq.0) then
    write(150,'("  interstitial time (seconds) : ",F12.2)')timeistl
    write(150,'("    muffin-tin time (seconds) : ",F12.2)')timemt
    call flushifc(150)
  endif
  
  write(160)zrhofc(1:ngvecme,1:num_nnp(ik),1)
enddo !ik
close(160)

deallocate(zrhofc)

#endif

deallocate(ikptiprocnr)
deallocate(ikq)
deallocate(vgq0c)
deallocate(gq0)
deallocate(tpgq0)
deallocate(sfacgq0)
deallocate(ylmgq0) 
deallocate(num_nnp)
deallocate(nnp)
deallocate(docc)

if (iproc.eq.0) then
  write(150,*)
  write(150,'("Done.")')
  call flushifc(150)
endif

return
end


subroutine calc_uuj(uuj,lmaxexp,gq0,ngvecme)
use modmain
implicit none
! arguments
integer, intent(in) :: lmaxexp
integer, intent(in) :: ngvecme
real(8), intent(in) :: gq0(ngvecme)
real(8), intent(out) :: uuj(0:lmaxvr,0:lmaxvr,0:lmaxexp,nrfmax,nrfmax,natmtot,ngvecme)
! local variables
integer ia,is,ias,l,io,ilo,ig,l1,l2,l3,io1,io2,ir
real(8), allocatable :: jlgq0r(:,:,:,:)
integer ordl(0:lmaxvr)
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

subroutine zrhoftit(ngvecme,gvecme1,max_num_nnp,num_nnp,nnp,ngknr1,ngknr2,igkignr1,igkignr2, &
  igkq,wfsvit1,wfsvit2,zrhofc)
use modmain
implicit none
! arguments
integer, intent(in) :: ngvecme
integer, intent(in) :: gvecme1
integer, intent(in) :: max_num_nnp
integer, intent(in) :: num_nnp
integer, intent(in) :: nnp(max_num_nnp,3)
integer, intent(in) :: ngknr1
integer, intent(in) :: ngknr2
integer, intent(in) :: igkq
integer, intent(in) :: igkignr1(ngkmax)
integer, intent(in) :: igkignr2(ngkmax)
complex(8), intent(in) :: wfsvit1(nmatmax,nstsv,nspinor)
complex(8), intent(in) :: wfsvit2(nmatmax,nstsv,nspinor)
complex(8), intent(inout) :: zrhofc(ngvecme,max_num_nnp)

complex(8), allocatable :: mit(:,:)
complex(8), allocatable :: a(:,:,:) 

integer is,ia,ias,ig,ig1,ig2,ist1,ist2,i,i1,i2,ispn,ispn2,istfv,nst1
integer iv3g(3)
real(8) v1(3),v2(3),tp3g(2),len3g
complex(8) sfac3g(natmtot)
complex(8) zt1

allocate(a(ngknr2,nstsv,nspinor))

! compute coefficients for required spin(s) only
!  for both spins compute nstfv+nstfv=nstsv states
!  for one (up or down spin) compute nstfv states
!if (spin_me.ne.3) then
!  nst1=nstfv
!else
!  nst1=nstsv
!endif
! for first spin or for two spins start summation from 1
!  for second spin start summation from nstfv+1 (skip spin up block)
!if (spin_me.eq.1.or.spin_me.eq.3) then
!  ispn=1
!else
!  ispn=2
!endif

allocate(mit(ngknr1,ngknr2))

do ig=1,ngvecme
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
    if (lrtype.eq.1.or.lrtype.eq.2) then
      ispn2=3-ispn
    endif
    do i=1,num_nnp
      ist1=nnp(i,1)
      ist2=nnp(i,2)
      do ig2=1,ngknr2
        zrhofc(ig,i)=zrhofc(ig,i)+wfsvit2(ig2,ist2,ispn2)*a(ig2,ist1,ispn)
      enddo
    enddo
  enddo
enddo
deallocate(mit,a)
return
end
        
subroutine zrhoftmt(ngvecme,max_num_nnp,num_nnp,nnp, &
  wfsvmt1,wfsvmt2,ngumax,ngu,gu,igu,zrhofc)
use modmain
implicit none
! arguments
integer, intent(in) :: ngvecme
integer, intent(in) :: max_num_nnp
integer, intent(in) :: num_nnp
integer, intent(in) :: nnp(max_num_nnp,3)
integer, intent(in) :: ngumax
integer, intent(in) :: ngu(natmtot,ngvecme)
integer, intent(in) :: igu(4,ngumax,natmtot,ngvecme)
complex(4), intent(in) :: gu(ngumax,natmtot,ngvecme)
complex(8), intent(in) :: wfsvmt1(lmmaxvr,nrfmax,natmtot,nstsv,nspinor)
complex(8), intent(in) :: wfsvmt2(lmmaxvr,nrfmax,natmtot,nstsv,nspinor)
complex(8), intent(inout) :: zrhofc(ngvecme,max_num_nnp)
! local variables
integer ig,i,j,ist1,ist2,ias,io1,io2,lm1,lm2,ispn,ispn2 
complex(8) a1(lmmaxvr,nrfmax),a2(lmmaxvr,nrfmax)

do ig=1,ngvecme
  do ispn=1,nspinor
    if (lrtype.eq.0) then
      ispn2=ispn
    endif
    if (lrtype.eq.1.or.lrtype.eq.2) then
      ispn2=3-ispn
    endif
  do i=1,num_nnp
    ist1=nnp(i,1)
    ist2=nnp(i,2)
    do ias=1,natmtot
      a1=dconjg(wfsvmt1(:,:,ias,ist1,ispn))
      a2=wfsvmt2(:,:,ias,ist2,ispn2)
      do j=1,ngu(ias,ig)
        lm1=igu(1,j,ias,ig)
        lm2=igu(2,j,ias,ig)
        io1=igu(3,j,ias,ig)
        io2=igu(4,j,ias,ig)
        zrhofc(ig,i)=zrhofc(ig,i)+a1(lm1,io1)*a2(lm2,io2)*gu(j,ias,ig)
      enddo
    enddo !ias 
  enddo !i
  enddo
enddo !ig    

return
end


subroutine getnnp(req,ikq,occsvnr,max_num_nnp,num_nnp,nnp,docc)
use modmain
implicit none
! arguments
logical, intent(in) :: req
integer, intent(in) :: ikq(nkptnr,2)
real(8), intent(in) :: occsvnr(nstsv,nkptnr)
integer, intent(inout) :: max_num_nnp
integer, intent(out) :: num_nnp(nkptnr)
integer, intent(out) :: nnp(max_num_nnp,3,nkptnr)
real(8), intent(out) :: docc(max_num_nnp,nkptnr)

integer band1,band2
integer i,ik,jk,istfv1,istfv2,ispn1,ispn2,ist1,ist2,ikloc
logical laddme,ldocc
real(8) d1
integer, external :: iknrglob

if (bndme1.eq.-1) then
  band1=1
  band2=nstfv
else
  band1=bndme1
  band2=bndme2
endif
if (req.and.iproc.eq.0) then
  write(150,*)
  write(150,'("Band interval: ",2I4)')band1,band2
endif
if (req) max_num_nnp=0
do ikloc=1,nkptnrloc(iproc)
  ik=iknrglob(ikloc)
  jk=ikq(ik,1)
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
            nnp(i,1,ikloc)=ist1
            nnp(i,2,ikloc)=ist2
            nnp(i,3,ikloc)=ispn1
            docc(i,ikloc)=d1
          endif
        endif
      enddo !istfv2
      enddo !istfv1
    enddo !ispn2
  enddo !ispn1
  if (.not.req) num_nnp(ikloc)=i
  if (req) max_num_nnp=max(max_num_nnp,i)
enddo !ikloc

return
end

subroutine getgu(req,ngvecme,lmaxexp,uuj,ylmgq0,sfacgq0,ngumax,ngu,gu,igu)
use modmain
implicit none
! arguments
logical, intent(in) :: req
integer, intent(in) :: ngvecme
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


