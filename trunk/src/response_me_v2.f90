subroutine response_me_v2(ivq0m,ngvec_me)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
! q-vector in k-mesh coordinates
integer, intent(in) :: ivq0m(3)
! number of G-vectors for matrix elements calculation
integer, intent(in) :: ngvec_me
 
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
!              2: index of equivalent k-point in irreducible part of BZ
!              3: index of equivalent k'-point in irreducible part of BZ
!              3: index of K-vector which brings k+q to first BZ
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
! Bessel functions j_l(|G+q|x)
real(8), allocatable :: jlgq0r(:,:,:,:)

! allocatable arrays
integer, allocatable :: ngknr(:)
integer, allocatable :: igkignr(:,:)
integer, allocatable :: igkignr2(:)
real(8), allocatable :: vgklnr(:,:,:)
real(8), allocatable :: vgklnr2(:,:)
real(8), allocatable :: vgkcnr(:,:)
real(8), allocatable :: gknr(:,:)
real(8), allocatable :: gknr2(:)
real(8), allocatable :: tpgknr(:,:,:)
real(8), allocatable :: tpgknr2(:,:)
real(8), allocatable :: occsvnr(:,:)
real(8), allocatable :: docc(:,:)
complex(8), allocatable :: sfacgknr(:,:,:)
complex(8), allocatable :: sfacgknr2(:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhofc(:,:,:)
complex(8), allocatable :: zrhofc0(:)
complex(8), allocatable :: evecfv2(:,:,:)
complex(8), allocatable :: evecsv2(:,:)

complex(8), allocatable :: apwcoeffloc(:,:,:,:,:,:)
complex(8), allocatable :: locoeffloc(:,:,:,:,:,:)
complex(8), allocatable :: apwcoeff2(:,:,:,:,:)
complex(8), allocatable :: locoeff2(:,:,:,:,:)

integer i,j,i1,ik,jk,ig,is,ir,ikstep,ist1,ist2,ispn,ikloc
integer ngknr2
real(8) vkq0l(3),t1,jl(0:lmaxvr)
integer ivg1(3)
real(8) cpu0,cpu1

! for parallel execution
integer, allocatable :: nkptlocnr(:)
integer, allocatable :: ikptlocnr(:,:)
integer, allocatable :: ikptiproc(:)
integer, allocatable :: ikptiprocnr(:)
integer, allocatable :: isend(:,:,:)
integer tag,req,ierr
integer, allocatable :: status(:)
character*100 :: fname

! external functions
real(8), external :: r3taxi
complex(8), external :: zfint

! read the density and potentials from file
call readstate

! read Fermi energy from file
call readfermi

! find the new linearisation energies
call linengy

! generate the APW radial functions
call genapwfr

! generate the local-orbital radial functions
call genlofr

if (iproc.eq.0) then
  write(150,'("Calculation of matrix elements <n,k|e^{-i(G+q)x}|n'',k+q>")')
  if (spinpol) then
    write(150,*)
    write(150,'("Spin-polarized calculation")')
    if (spin_me.eq.1) write(150,'(" calculation of matrix elements for spin up")')
    if (spin_me.eq.2) write(150,'(" calculation of matrix elements for spin dn")')
    if (spin_me.eq.3) write(150,'(" calculation of matrix elements for both spins")')
  endif
  write(150,*)
  write(150,'("Number of G-shells  : ",I4)')ngsh_me
  write(150,'("Number of G-vectors : ",I4)')ngvec_me
endif

allocate(nkptlocnr(0:nproc-1))
allocate(ikptlocnr(0:nproc-1,2))
call splitk(nkptnr,nproc,nkptlocnr,ikptlocnr)
allocate(ikptiproc(nkpt))
allocate(ikptiprocnr(nkptnr))
do i=0,nproc-1
  ikptiproc(ikptloc(i,1):ikptloc(i,2))=i
  ikptiprocnr(ikptlocnr(i,1):ikptlocnr(i,2))=i
enddo
!if (iproc.eq.0.and.nproc.gt.1) then
!  write(150,*)
!  write(150,'("Reduced k-points division:")')
!  write(150,'(" iproc  first k   last k   nkpt ")')
!  write(150,'(" ------------------------------ ")')
!  do i=0,nproc-1
!    write(150,'(1X,I4,4X,I4,5X,I4,5X,I4)')i,ikptloc(i,1),ikptloc(i,2),nkptloc(i)
!  enddo
!  write(150,*)
!  write(150,'("Non-reduced k-points division:")')
!   write(150,'(" iproc  first k   last k   nkpt ")')
!  write(150,'(" ------------------------------ ")')
!  do i=0,nproc-1
!    write(150,'(1X,I4,4X,I4,5X,I4,5X,I4)')i,ikptlocnr(i,1),ikptlocnr(i,2),nkptlocnr(i)
!  enddo
!endif

! generate G+k vectors for entire BZ (this is required to compute 
!   wave-functions at each k-point)
allocate(vgklnr(3,ngkmax,nkptlocnr(iproc)))
allocate(vgkcnr(3,ngkmax))
allocate(gknr(ngkmax,nkptlocnr(iproc)))
allocate(tpgknr(2,ngkmax,nkptlocnr(iproc)))
allocate(ngknr(nkptlocnr(iproc)))
allocate(sfacgknr(ngkmax,natmtot,nkptlocnr(iproc)))
allocate(igkignr(ngkmax,nkptlocnr(iproc)))
do ikloc=1,nkptlocnr(iproc)
  ik=ikptlocnr(iproc,1)+ikloc-1
  call gengpvec(vklnr(1,ik),vkcnr(1,ik),ngknr(ikloc),igkignr(1,ikloc), &
    vgklnr(1,1,ikloc),vgkcnr,gknr(1,ikloc),tpgknr(1,1,ikloc))
  call gensfacgp(ngknr(ikloc),vgkcnr,ngkmax,sfacgknr(1,1,ikloc))
enddo

allocate(ikq(nkptnr,4))
allocate(vgq0c(3,ngvec))
allocate(gq0(ngvec))
allocate(tpgq0(2,ngvec))
allocate(sfacgq0(ngvec,natmtot))
allocate(ylmgq0(lmmaxvr,ngvec)) 
allocate(jlgq0r(nrcmtmax,0:lmaxvr,nspecies,ngvec))
allocate(occsvnr(nstsv,nkptnr))

! q-vector in lattice coordinates
do i=1,3
  vq0l(i)=1.d0*ivq0m(i)/ngridk(i)
enddo

! find G-vector which brings q0 to first BZ
vgq0l(:)=floor(vq0l(:))

! reduce q0 vector to first BZ
vq0rl(:)=vq0l(:)-vgq0l(:)

! check if we have enough G-shells to bring q-vector back to first BZ
do ig=1,ngvec_me
  if (sum(abs(vgq0l(:)-ivg(:,ig))).eq.0) then
    igq0=ig
    goto 20
  endif
enddo
write(*,*)
write(*,'("Bug(response_me): not enough G-vectors to reduce q-vector &
  &to first BZ")')
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
  write(150,'("reduced q-vector (lat.coord.)                : ",&
    & 3G18.10)')vq0rl
  write(150,'("reduced q-vector (Cart.coord.) [a.u.]        : ",&
    & 3G18.10)')vq0rc
endif

! find k+q and reduce them to first BZ (this is required to utilize the 
!   periodical property of Bloch-states: |k>=|k+K>, where K is any vector 
!   of the reciprocal lattice)
!if (iproc.eq.0) then
!  write(150,*)
!  write(150,'(3X,"ik",10X,"k",19X,"k+q",16X,"K",12X,"k''=k+q-K",8X,"jk")')
!  write(150,'(85("-"))')
!endif
do ik=1,nkptnr
! k+q vector
  vkq0l(:)=vklnr(:,ik)+vq0rl(:)+1d-10
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
  ikq(ik,4)=ivgig(ivg1(1),ivg1(2),ivg1(3))
!  if (iproc.eq.0) then
!    write(150,'(I4,2X,3F6.2,2X,3F6.2,2X,3I4,2X,3F6.2,2X,I4)') &
!    ik,vklnr(:,ik),vkq0l+ivg1,ivg1,vkq0l,ikq(ik,1)
!  endif
! for parallel execution
  if (ismpi) then
    call findkpt(vklnr(1,ik),i,ikq(ik,2))
    call findkpt(vklnr(1,ikq(ik,1)),i,ikq(ik,3))
  endif
enddo

! get occupancy of states
if (iproc.eq.0) then 
  do ik=1,nkptnr
    call getoccsv(vklnr(1,ik),occsvnr(1,ik))
  enddo
endif
#ifdef _MPI_
call mpi_bcast(occsvnr,nstsv*nkptnr,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif

! setup n,n' stuff
! first, find the maximum size of nnp array
max_num_nnp=0
allocate(num_nnp(nkptnr))
do ik=1,nkptnr
  jk=ikq(ik,1)
  i1=0
  do ispn=1,nspinor
    do i=1,nstfv
      ist1=i+(ispn-1)*nstfv
      do j=1,nstfv
        ist2=j+(ispn-1)*nstfv
        if ((ispn.eq.spin_me.or.spin_me.eq.3) .and. &
	    abs(occsvnr(ist1,ik)-occsvnr(ist2,jk)).gt.1d-10) i1=i1+1
      enddo
    enddo
  enddo
  num_nnp(ik)=i1
  max_num_nnp=max(max_num_nnp,i1)
enddo
if (iproc.eq.0) then
  write(150,*)
  write(150,'("Maximum number of interband transitions: ",I5)')max_num_nnp
endif
allocate(nnp(max_num_nnp,3,nkptnr))
allocate(docc(nkptnr,max_num_nnp))
! second, setup the nnp array
do ik=1,nkptnr
  jk=ikq(ik,1)
  i1=0
  do ispn=1,nspinor
    do i=1,nstfv
      ist1=i+(ispn-1)*nstfv
      do j=1,nstfv
        ist2=j+(ispn-1)*nstfv
	if ((ispn.eq.spin_me.or.spin_me.eq.3) .and. &
	    abs(occsvnr(ist1,ik)-occsvnr(ist2,jk)).gt.1d-10) then
          i1=i1+1
          nnp(i1,1,ik)=ist1
          nnp(i1,2,ik)=ist2
	  if (spin_me.eq.3) then
	    nnp(i1,3,ik)=ispn
	  else
            nnp(i1,3,ik)=1
	  endif
          docc(ik,i1)=occsvnr(ist1,ik)-occsvnr(ist2,jk)
        endif
      enddo
    enddo
  enddo
enddo

! generate G+q' vectors, where q' is reduced q-vector
do ig=1,ngvec_me
  vgq0c(:,ig)=vgc(:,ig)+vq0rc(:)
! get spherical coordinates and length of G+q'
  call sphcrd(vgq0c(:,ig),gq0(ig),tpgq0(:,ig))
! generate spherical harmonics for G+q'
  call genylm(lmaxvr,tpgq0(:,ig),ylmgq0(:,ig))
enddo

! generate structure factor for G+q' vectors
call gensfacgp(ngvec_me,vgq0c,ngvec,sfacgq0)
  
! generate Bessel functions j_l(|G+q'|x)
do ig=1,ngvec_me
  do is=1,nspecies
    do ir=1,nrcmt(is)
      t1=gq0(ig)*rcmt(ir,is)
      call sbessel(lmaxvr,t1,jl)
      jlgq0r(ir,:,is,ig)=jl(:)
    enddo
  enddo
enddo

write(fname,'("ZRHOFC[",I4.3,",",I4.3,",",I4.3,"].OUT")') &
  ivq0m(1),ivq0m(2),ivq0m(3)

if (iproc.eq.0) then
  open(160,file=trim(fname),form='unformatted',status='replace')
  write(160)nkptnr,ngsh_me,ngvec_me,max_num_nnp,igq0
  write(160)nspinor,spin_me
  write(160)vq0l(1:3)
  write(160)vq0rl(1:3)
  write(160)vq0c(1:3)
  write(160)vq0rc(1:3)
  do ik=1,nkptnr
    write(160)ikq(ik,1)
    write(160)num_nnp(ik)
    write(160)nnp(1:num_nnp(ik),1:3,ik)
    write(160)docc(ik,1:num_nnp(ik))
  enddo
  close(160)
endif

allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
allocate(zrhofc0(ngvec_me))
  
! different implementation for parallel and serial execution
#ifdef _MPI_

allocate(evecfvloc(nmatmax,nstfv,nspnfv,nkptlocnr(iproc)))
allocate(evecsvloc(nstsv,nstsv,nkptlocnr(iproc)))
allocate(apwcoeffloc(nstsv,nspinor,apwordmax,lmmaxvr,natmtot,nkptlocnr(iproc)))
allocate(locoeffloc(nstsv,nspinor,nlomax,lmmaxvr,natmtot,nkptlocnr(iproc)))
allocate(status(MPI_STATUS_SIZE))
allocate(apwcoeff2(nstsv,nspinor,apwordmax,lmmaxvr,natmtot))
allocate(locoeff2(nstsv,nspinor,nlomax,lmmaxvr,natmtot))
allocate(evecfv2(nmatmax,nstfv,nspnfv))
allocate(evecsv2(nstsv,nstsv))
allocate(igkignr2(ngkmax))
allocate(zrhofc(ngvec_me,max_num_nnp,nkptlocnr(iproc)))

! read all eigen-vectors to the memory
if (iproc.eq.0) then
  write(150,*)
  write(150,'("Reading eigen-vectors")')
  call flushifc(150)
endif
do i=0,nproc-1
  if (iproc.eq.i) then
    do ikloc=1,nkptlocnr(iproc)
      ik=ikptlocnr(iproc,1)+ikloc-1
      call getevecfv(vklnr(1,ik),vgklnr(1,1,ikloc),evecfvloc(1,1,1,ikloc))
      call getevecsv(vklnr(1,ik),evecsvloc(1,1,ikloc))
    enddo 
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierr)
enddo
if (iproc.eq.0) then
  write(150,'("Done.")')
  call flushifc(150)
endif

! calculate muffin-tin coefficients
if (iproc.eq.0) then
  write(150,*)
  write(150,'("Calculating muffin-tin coefficients")')
  call flushifc(150)
endif
do ikloc=1,nkptlocnr(iproc)
  call match(ngknr(ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc),sfacgknr(1,1,ikloc),apwalm)
  call tcoeff(ngknr(ikloc),apwalm,evecfvloc(1,1,1,ikloc),evecsvloc(1,1,ikloc),&
    apwcoeffloc(1,1,1,1,1,ikloc),locoeffloc(1,1,1,1,1,ikloc))
enddo
if (iproc.eq.0) then
  write(150,'("Done.")')
  call flushifc(150)
endif

! find indexes of k-points to send and receive
allocate(isend(nkptlocnr(0),0:nproc-1,2))
isend=-1
do ikstep=1,nkptlocnr(0)
  do i=0,nproc-1
    ik=ikptlocnr(i,1)+ikstep-1
    if (ikstep.le.nkptlocnr(i)) then
! for non reduced k' point find the proc and local index
      isend(ikstep,i,1)=ikptiprocnr(ikq(ik,1))
      isend(ikstep,i,2)=ikq(ik,1)-ikptlocnr(isend(ikstep,i,1),1)+1
    endif
  enddo
enddo

if (iproc.eq.0) then
  write(150,*)
  write(150,'("Starting k-point loop")')
endif
do ikstep=1,nkptlocnr(0)
  if (iproc.eq.0) then
    write(150,'("k-step ",I4," out of ",I4)')ikstep,nkptlocnr(0)
    call flushifc(150)
    call cpu_time(cpu0)
  endif
! find the i-th proc to which the current iproc should send data
  do i=0,nproc-1
    if (isend(ikstep,i,1).eq.iproc.and.iproc.ne.i) then
      tag=(ikstep*nproc+i)*10
      ik=isend(ikstep,i,2)
      call mpi_isend(evecfvloc(1,1,1,ik),                    &
        nmatmax*nstfv*nspnfv,                                &
        MPI_DOUBLE_COMPLEX,i,tag,MPI_COMM_WORLD,req,ierr)
      tag=tag+1
      call mpi_isend(evecsvloc(1,1,ik),                      &
        nstsv*nstsv,                                         &
        MPI_DOUBLE_COMPLEX,i,tag,MPI_COMM_WORLD,req,ierr)
      tag=tag+1
      call mpi_isend(apwcoeffloc(1,1,1,1,1,ik),              &
        nstsv*nspinor*apwordmax*lmmaxvr*natmtot,             &
	MPI_DOUBLE_COMPLEX,i,tag,MPI_COMM_WORLD,req,ierr)
      tag=tag+1
      call mpi_isend(locoeffloc(1,1,1,1,1,ik),               &
        nstsv*nspinor*nlomax*lmmaxvr*natmtot,                & 
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
      call mpi_recv(evecfv2,nmatmax*nstfv*nspnfv,MPI_DOUBLE_COMPLEX, &
        isend(ikstep,iproc,1),tag,MPI_COMM_WORLD,status,ierr)
      tag=tag+1
      call mpi_recv(evecsv2,nstsv*nstsv,MPI_DOUBLE_COMPLEX,          &
        isend(ikstep,iproc,1),tag,MPI_COMM_WORLD,status,ierr)
      tag=tag+1
      call mpi_recv(apwcoeff2,nstsv*nspinor*apwordmax*lmmaxvr*natmtot, &
        MPI_DOUBLE_COMPLEX,isend(ikstep,iproc,1),tag,MPI_COMM_WORLD,status,ierr)
      tag=tag+1
      call mpi_recv(locoeff2,nstsv*nspinor*nlomax*lmmaxvr*natmtot, &
        MPI_DOUBLE_COMPLEX,isend(ikstep,iproc,1),tag,MPI_COMM_WORLD,status,ierr)
      tag=tag+1
      call mpi_recv(ngknr2,1,MPI_INTEGER,isend(ikstep,iproc,1),tag,MPI_COMM_WORLD,status,ierr)
      tag=tag+1
      call mpi_recv(igkignr2,ngkmax,MPI_INTEGER,isend(ikstep,iproc,1),tag,MPI_COMM_WORLD,status,ierr)
    else
      ik=isend(ikstep,iproc,2)
      evecfv2(:,:,:)=evecfvloc(:,:,:,ik)
      evecsv2(:,:)=evecsvloc(:,:,ik)
      apwcoeff2(:,:,:,:,:)=apwcoeffloc(:,:,:,:,:,ik)
      locoeff2(:,:,:,:,:)=locoeffloc(:,:,:,:,:,ik)
      ngknr2=ngknr(ik)
      igkignr2(:)=igkignr(:,ik)
    endif
  endif

  call mpi_barrier(MPI_COMM_WORLD,ierr)
  if (iproc.eq.0) then
    write(150,'("  OK send and recieve")')
    call flushifc(150)
  endif
  
  if (ikstep.le.nkptlocnr(iproc)) then
    ik=ikptlocnr(iproc,1)+ikstep-1
    jk=ikq(ik,1)
    
! calculate interstitial contribution for all combinations of n,n'
    call zrhoftistl(ngvec_me,max_num_nnp,num_nnp(ik),nnp(1,1,ik),ngknr(ikstep),ngknr2, &
      igkignr(1,ikstep),igkignr2,ikq(ik,4),evecfvloc(1,1,1,ikstep),evecsvloc(1,1,ikstep), &
      evecfv2,evecsv2,zrhofc(1,1,ikstep))
    
    if (iproc.eq.0) then
      write(150,'("  Intertitial contribution is done")')
      call flushifc(150)
    endif

    do i=1,num_nnp(ik)
      ist1=nnp(i,1,ik)
      ist2=nnp(i,2,ik)
      call zrhomt0(ist1,ist2,apwcoeffloc(1,1,1,1,1,ikstep),locoeffloc(1,1,1,1,1,ikstep),&
        apwcoeff2,locoeff2,zrhomt)
      call zrhoftmt(zrhomt,jlgq0r,ylmgq0,sfacgq0,ngvec_me,zrhofc0)
      zrhofc(:,i,ikstep)=zrhofc(:,i,ikstep)+zrhofc0(:)
    enddo
  endif ! (ikstep.le.nkptlocnr(iproc))
  
  if (iproc.eq.0) then
    call cpu_time(cpu1)
    write(150,'("  step time (seconds) : ",F12.2)')cpu1-cpu0
    call flushifc(150)
  endif

enddo !ikstep

call mpi_barrier(MPI_COMM_WORLD,ierr)

do i=0,nproc-1
  if (i.eq.iproc) then
    open(160,file=trim(fname),form='unformatted',status='old',position='append')
    do ikstep=1,nkptlocnr(iproc)
      ik=ikptlocnr(iproc,1)+ikstep-1
      write(160)ik,ikq(ik,1)
      write(160)zrhofc(1:ngvec_me,1:num_nnp(ik),ikstep)
    enddo !ikstep
    close(160)
  endif !i.eq.iproc
  call mpi_barrier(MPI_COMM_WORLD,ierr)
enddo 

deallocate(evecfvloc)
deallocate(evecsvloc)
deallocate(apwcoeffloc)
deallocate(locoeffloc)
deallocate(status)
deallocate(apwcoeff2)
deallocate(locoeff2)
deallocate(evecfv2)
deallocate(evecsv2)
deallocate(igkignr2)
deallocate(zrhofc)
deallocate(isend)
  
#else

allocate(evecfv1(nmatmax,nstfv,nspnfv))
allocate(evecsv1(nstsv,nstsv))
allocate(evecfv2(nmatmax,nstfv,nspnfv))
allocate(evecsv2(nstsv,nstsv))
allocate(zrhofc(ngvec_me,max_num_nnp,1))

allocate(apwcoeff1(nstsv,nspinor,apwordmax,lmmaxvr,natmtot))
allocate(apwcoeff2(nstsv,nspinor,apwordmax,lmmaxvr,natmtot))
allocate(locoeff1(nstsv,nspinor,nlomax,lmmaxvr,natmtot))
allocate(locoeff2(nstsv,nspinor,nlomax,lmmaxvr,natmtot))


open(160,file=trim(fname),form='unformatted',status='old',position='append')

write(150,*)
do ik=1,nkptnr
  write(150,'("k-point ",I4," out of ",I4)')ik,nkptnr
  call flushifc(150)
  
  jk=ikq(ik,1)
  
  write(160)ik,jk
  
! generate data at k
  call getevecfv(vklnr(1,ik),vgklnr(1,1,ik),evecfv1)
  call getevecsv(vklnr(1,ik),evecsv1) 
  call match(ngknr(ik),gknr(1,ik),tpgknr(1,1,ik),sfacgknr(1,1,ik),apwalm)
  call tcoeff(ngknr(ik),apwalm,evecfv1,evecsv1,apwcoeff1,locoeff1)

! generate data at k'=k+q-K
  call getevecfv(vklnr(1,jk),vgklnr(1,1,jk),evecfv2)
  call getevecsv(vklnr(1,jk),evecsv2) 
  call match(ngknr(jk),gknr(1,jk),tpgknr(1,1,jk),sfacgknr(1,1,jk),apwalm)
  call tcoeff(ngknr(jk),apwalm,evecfv2,evecsv2,apwcoeff2,locoeff2)

! calculate interstitial contribution for all combinations of n,n'
  call zrhoftistl(ngvec_me,max_num_nnp,num_nnp(ik),nnp(1,1,ik),ngknr(ik),ngknr(jk), &
    igkignr(1,ik),igkignr(1,jk),ikq(ik,4),evecfv1,evecsv1,evecfv2,evecsv2,zrhofc)
  
  do i=1,num_nnp(ik)
    ist1=nnp(i,1,ik)
    ist2=nnp(i,2,ik)
    call zrhomt0(ist1,ist2,apwcoeff1,locoeff1,apwcoeff2,locoeff2,zrhomt)
    call zrhoftmt(zrhomt,jlgq0r,ylmgq0,sfacgq0,ngvec_me,zrhofc1(1,1))
    zrhofc(:,i,1)=zrhofc(:,i,1)+zrhofc1(:,1)
  enddo
  
  write(160)zrhofc(1:ngvec_me,1:num_nnp(ik),1)
enddo !ik
close(160)

deallocate(evecfv1)
deallocate(evecsv1)
deallocate(evecfv2)
deallocate(evecsv2)
deallocate(zrhofc)
deallocate(apwcoeff1)
deallocate(apwcoeff2)
deallocate(locoeff1)
deallocate(locoeff2)

#endif

deallocate(nkptlocnr)
deallocate(ikptlocnr)
deallocate(ikptiproc)
deallocate(ikptiprocnr)
deallocate(vgklnr)
deallocate(vgkcnr)
deallocate(gknr)
deallocate(tpgknr)
deallocate(ngknr)
deallocate(sfacgknr)
deallocate(igkignr)
deallocate(apwalm)
deallocate(zrhomt)
deallocate(ikq)
deallocate(vgq0c)
deallocate(gq0)
deallocate(tpgq0)
deallocate(sfacgq0)
deallocate(ylmgq0) 
deallocate(jlgq0r)
deallocate(occsvnr)
deallocate(zrhofc0)
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


subroutine tcoeff(ngp,apwalm,evecfv,evecsv,apwcoeff,locoeff)
use modmain
implicit none
! arguments
integer, intent(in) :: ngp
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
complex(8), intent(out) :: apwcoeff(nstsv,nspinor,apwordmax,lmmaxvr,natmtot)
complex(8), intent(out) :: locoeff(nstsv,nspinor,nlomax,lmmaxvr,natmtot)
! local variables
integer i,j,l,m,ispn,istfv,is,ia,ias,lm,ig,i1,io,ilo

apwcoeff=dcmplx(0.d0,0.d0)
locoeff=dcmplx(0.d0,0.d0)
do j=1,nstsv
  i=0
  do ispn=1,nspinor
    do istfv=1,nstfv
      i=i+1
      do is=1,nspecies
        do ia=1,natoms(is)
          ias=idxas(ia,is)
! apw coefficients
          do l=0,lmaxvr
            do m=-l,l
              lm=idxlm(l,m)
              do ig=1,ngp
                do io=1,apword(l,is)
                  apwcoeff(j,ispn,io,lm,ias)=apwcoeff(j,ispn,io,lm,ias) + &
                    evecsv(i,j)*evecfv(ig,istfv)*apwalm(ig,io,lm,ias)
                enddo !io
              enddo !ig
            enddo !m
          enddo !l
! local orbital coefficients     
          do ilo=1,nlorb(is)
            l=lorbl(ilo,is)
            if (l.le.lmaxvr) then
              do m=-l,l
                lm=idxlm(l,m)
                i1=ngp+idxlo(lm,ilo,ias)
                locoeff(j,ispn,ilo,lm,ias)=locoeff(j,ispn,ilo,lm,ias) + &
                  evecsv(i,j)*evecfv(i1,istfv)
              enddo !m
            endif
          enddo !ilo
        enddo !ia
      enddo !is
    enddo !ist
  enddo !ispn
enddo !j

return
end

subroutine zrhomt0(ist1,ist2,apwcoeff1,locoeff1,apwcoeff2,locoeff2,zrhomt)
use modmain
implicit none
! arguments
integer, intent(in) :: ist1
integer, intent(in) :: ist2
complex(8), intent(in) :: apwcoeff1(nstsv,nspinor,apwordmax,lmmaxvr,natmtot)
complex(8), intent(in) :: locoeff1(nstsv,nspinor,nlomax,lmmaxvr,natmtot)
complex(8), intent(in) :: apwcoeff2(nstsv,nspinor,apwordmax,lmmaxvr,natmtot)
complex(8), intent(in) :: locoeff2(nstsv,nspinor,nlomax,lmmaxvr,natmtot)
complex(8), intent(out) :: zrhomt(lmmaxvr,nrcmtmax,natmtot)
! local variables
integer ispn,is,ia,ias,l,m,lm,io,ilo,ir,irc
integer l1,m1,lm1,l2,m2,lm2,l3,m3,lm3
! automatic arrays
complex(8) wfmt1(lmmaxvr,nrcmtmax)
complex(8) wfmt2(lmmaxvr,nrcmtmax)
complex(8) wfmt1t(lmmaxvr,nrcmtmax)
complex(8) wfmt2t(lmmaxvr,nrcmtmax)
complex(8) zrhot(lmmaxvr,nrcmtmax)

do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    zrhot=dcmplx(0.d0,0.d0)
    do ispn=1,nspinor
      wfmt1=dcmplx(0.d0,0.d0)
      wfmt2=dcmplx(0.d0,0.d0)
! add apw contribution
      do l=0,lmaxvr       
        do m=-l,l
          lm=idxlm(l,m)
          do io=1,apword(l,is)
            if (abs(apwcoeff1(ist1,ispn,io,lm,ias)).gt.1d-12) then
              irc=0
              do ir=1,nrmt(is),lradstp
                irc=irc+1
                wfmt1(lm,irc)=wfmt1(lm,irc)+apwcoeff1(ist1,ispn,io,lm,ias) * &
                  apwfr(ir,1,io,l,ias)
              enddo !ir
            endif
            if (abs(apwcoeff2(ist2,ispn,io,lm,ias)).gt.1d-12) then
              irc=0
              do ir=1,nrmt(is),lradstp
                irc=irc+1
                wfmt2(lm,irc)=wfmt2(lm,irc)+apwcoeff2(ist2,ispn,io,lm,ias) * &
                  apwfr(ir,1,io,l,ias)
              enddo !ir
            endif
          enddo !io
        enddo !m
      enddo !l
! add local orbital contribution     
      do ilo=1,nlorb(is)
        l=lorbl(ilo,is)
        if (l.le.lmaxvr) then
          do m=-l,l
            lm=idxlm(l,m)
            if (abs(locoeff1(ist1,ispn,ilo,lm,ias)).gt.1d-12) then
              irc=0
              do ir=1,nrmt(is),lradstp
                irc=irc+1
                wfmt1(lm,irc)=wfmt1(lm,irc)+locoeff1(ist1,ispn,ilo,lm,ias) * &
                  lofr(ir,1,ilo,ias)
              enddo
            endif
            if (abs(locoeff2(ist2,ispn,ilo,lm,ias)).gt.1d-12) then
              irc=0
              do ir=1,nrmt(is),lradstp
                irc=irc+1
                wfmt2(lm,irc)=wfmt2(lm,irc)+locoeff2(ist2,ispn,ilo,lm,ias) * &
                  lofr(ir,1,ilo,ias)
              enddo
            endif
          enddo
        endif  
      enddo
      call zgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,zone,zbshtapw,lmmaxapw, &
           wfmt1,lmmaxvr,zzero,wfmt1t,lmmaxvr)
      call zgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,zone,zbshtapw,lmmaxapw, &
           wfmt2,lmmaxvr,zzero,wfmt2t,lmmaxvr)
      do irc=1,nrcmt(is)
        zrhot(:,irc)=zrhot(:,irc) + dconjg(wfmt1t(:,irc))*wfmt2t(:,irc)
      enddo
    enddo !ispn
    call zgemm('N','N',lmmaxvr,nrcmt(is),lmmaxvr,zone,zfshtvr,lmmaxvr,zrhot, &
      lmmaxvr,zzero,zrhomt(1,1,ias),lmmaxvr)
  enddo !ia
enddo !is

return
end

subroutine zrhoftistl(ngvec_me,max_num_nnp,num_nnp,nnp,ngknri,ngknrj,igkignri,igkignrj, &
  ig1,evecfvi,evecsvi,evecfvj,evecsvj,zrhofc)
use modmain
implicit none
! arguments
integer, intent(in) :: ngvec_me
integer, intent(in) :: max_num_nnp
integer, intent(in) :: num_nnp
integer, intent(in) :: nnp(max_num_nnp,3)
integer, intent(in) :: ngknri
integer, intent(in) :: ngknrj
integer, intent(in) :: ig1
integer, intent(in) :: igkignri(ngkmax)
integer, intent(in) :: igkignrj(ngkmax)
complex(8), intent(in) :: evecfvi(nmatmax,nstfv)
complex(8), intent(in) :: evecsvi(nstsv,nstsv)
complex(8), intent(in) :: evecfvj(nmatmax,nstfv)
complex(8), intent(in) :: evecsvj(nstsv,nstsv)
complex(8), intent(out) :: zrhofc(ngvec_me,max_num_nnp)


complex(8), allocatable :: mit(:,:)
complex(8), allocatable :: mit1(:,:)

integer is,ia,ias,ig,igi,igj,ist1,ist2,i,i1,i2,ispn
integer iv3g(3)
real(8) v1(3),v2(3),tp3g(2),len3g
complex(8) sfac3g(natmtot)
complex(8) zt1

allocate(mit(ngknri,ngknrj))
allocate(mit1(nstfv,nstfv))
zrhofc=dcmplx(0.d0,0.d0)
do ig=1,ngvec_me
  mit=dcmplx(0.d0,0.d0)
  do igi=1,ngknri
    do igj=1,ngknrj
      ! G1-G2+G+K
      iv3g(:)=ivg(:,igkignri(igi))-ivg(:,igkignrj(igj))+ivg(:,ig)+ivg(:,ig1)
      if (sum(abs(iv3g)).eq.0) mit(igi,igj)=dcmplx(1.d0,0.d0)
      v2(:)=1.d0*iv3g(:)
      call r3mv(bvec,v2,v1)
      call sphcrd(v1,len3g,tp3g)
      call gensfacgp(1,v1,1,sfac3g)
      do is=1,nspecies
        do ia=1,natoms(is)
	  ias=idxas(ia,is)
	  if (len3g.lt.1d-8) then
	    mit(igi,igj)=mit(igi,igj)-(fourpi/omega)*dconjg(sfac3g(ias))*(rmt(is)**3)/3.d0
	  else
	    mit(igi,igj)=mit(igi,igj)-(fourpi/omega)*dconjg(sfac3g(ias)) * &
	        (-(rmt(is)/len3g**2)*cos(len3g*rmt(is))+(1/len3g**3)*sin(len3g*rmt(is)))
	  endif
	enddo !ia
      enddo !is
    enddo
  enddo
  mit1=dcmplx(0.d0,0.d0)
  do ist1=1,nstfv
    do ist2=1,nstfv
      do igj=1,ngknrj
        zt1=dcmplx(0.d0,0.d0) 
        do igi=1,ngknri
          !mit1(ist1,ist2)=mit1(ist1,ist2)+dconjg(evecfvi(igi,ist1))*evecfvj(igj,ist2)*mit(igi,igj)
          zt1=zt1+dconjg(evecfvi(igi,ist1))*mit(igi,igj)
	enddo
	mit1(ist1,ist2)=mit1(ist1,ist2)+zt1*evecfvj(igj,ist2)
      enddo
    enddo
  enddo
      
  do i=1,num_nnp
    ist1=nnp(i,1)
    ist2=nnp(i,2)
    do ispn=1,nspinor
      do i2=1,nstfv
        zt1=dcmplx(0.d0,0.d0)
        do i1=1,nstfv
          !zrhofc(ig,i)=zrhofc(ig,i)+dconjg(evecsvi(i1+(ispn-1)*nstfv,ist1))*evecsvj(i2+(ispn-1)*nstfv,ist2)*mit1(i1,i2)
          zt1=zt1+dconjg(evecsvi(i1+(ispn-1)*nstfv,ist1))*mit1(i1,i2)
	enddo
	zrhofc(ig,i)=zrhofc(ig,i)+zt1*evecsvj(i2+(ispn-1)*nstfv,ist2)
      enddo
    enddo
!    zrhofc(ig,i)=mit1(ist1,ist2)
!    do igi=1,ngknri
!      do igj=1,ngknrj 
!        zrhofc(ig,i)=zrhofc(ig,i)+dconjg(evecfvi(igi,ist1))*evecfvj(igj,ist2)*mit(igi,igj)
!      enddo
!    enddo
  enddo
enddo
deallocate(mit,mit1)
return
end
