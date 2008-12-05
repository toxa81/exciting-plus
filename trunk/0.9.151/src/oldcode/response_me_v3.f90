subroutine response_me_v3(ivq0m,ngvec_me)
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
complex(8), allocatable :: evecfv1(:,:,:)
complex(8), allocatable :: evecsv1(:,:)
complex(8), allocatable :: evecfv2(:,:,:)
complex(8), allocatable :: evecsv2(:,:)

complex(8), allocatable :: acoeffloc(:,:,:,:,:,:)
complex(8), allocatable :: acoeff1(:,:,:,:,:)
complex(8), allocatable :: acoeff2(:,:,:,:,:)

integer i,j,i1,ik,jk,ig,is,ir,ikstep,ist1,ist2,ispn,ikloc,ilo,l,ias,io1,io2
integer ngknr2
real(8) vkq0l(3),t1,jl(0:lmaxvr)
integer ivg1(3)
real(8) cpu0,cpu1,timeistl,timemt
integer nlomaxl
integer mtord
integer, allocatable :: ltmp(:)
real(8), allocatable :: uuj(:,:,:,:,:,:)
complex(8), allocatable :: gu(:,:)
integer, allocatable :: igu(:,:,:)
integer, allocatable :: ngu(:)
integer ngumax
integer l1,m1,lm1,l2,m2,lm2,l3,m3,lm3

integer lmaxexp
integer lmmaxexp

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
real(8), external :: gaunt
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

lmaxexp=lmaxvr
lmmaxexp=(lmaxexp+1)**2

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
allocate(vgq0c(3,ngvec_me))
allocate(gq0(ngvec_me))
allocate(tpgq0(2,ngvec_me))
allocate(sfacgq0(ngvec_me,natmtot))
allocate(ylmgq0(lmmaxexp,ngvec_me)) 
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
allocate(num_nnp(nkptlocnr(iproc)))
max_num_nnp=nstfv**2
allocate(nnp(max_num_nnp,3,nkptlocnr(iproc)))
allocate(docc(max_num_nnp,nkptlocnr(iproc)))
do ikloc=1,nkptlocnr(iproc)
  ik=ikptlocnr(iproc,1)+ikloc-1
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
	  if (i1.gt.max_num_nnp) then
	    write(*,*)
	    write(*,'("Error(response_me): maximum number of nn'' combinations reached.")')
	    write(*,*)
	    call pstop
	  endif
          nnp(i1,1,ikloc)=ist1
          nnp(i1,2,ikloc)=ist2
	  if (spin_me.eq.3) then
	    nnp(i1,3,ikloc)=ispn
	  else
            nnp(i1,3,ikloc)=1
	  endif
          docc(i1,ikloc)=occsvnr(ist1,ik)-occsvnr(ist2,jk)
	endif
      enddo
    enddo
  enddo
  num_nnp(ikloc)=i1
enddo

if (iproc.eq.0) then
  write(150,*)
  write(150,'("Number of interband transitions in Gamma point : ",I5)')num_nnp(1)
endif

! generate G+q' vectors, where q' is reduced q-vector
do ig=1,ngvec_me
  vgq0c(:,ig)=vgc(:,ig)+vq0rc(:)
! get spherical coordinates and length of G+q'
  call sphcrd(vgq0c(:,ig),gq0(ig),tpgq0(:,ig))
! generate spherical harmonics for G+q'
  call genylm(lmaxexp,tpgq0(:,ig),ylmgq0(:,ig))
enddo

! generate structure factor for G+q' vectors
call gensfacgp(ngvec_me,vgq0c,ngvec_me,sfacgq0)

write(fname,'("ZRHOFC[",I4.3,",",I4.3,",",I4.3,"].OUT")') &
  ivq0m(1),ivq0m(2),ivq0m(3)

if (iproc.eq.0) then
  open(160,file=trim(fname),form='unformatted',status='replace')
  write(160)nkptnr,ngsh_me,ngvec_me,igq0
  write(160)nspinor,spin_me
  write(160)vq0l(1:3)
  write(160)vq0rl(1:3)
  write(160)vq0c(1:3)
  write(160)vq0rc(1:3)
  do ik=1,nkptnr
    write(160)ikq(ik,1)
!    write(160)num_nnp(ik)
!    write(160)nnp(1:num_nnp(ik),1:3,ik)
!    write(160)docc(ik,1:num_nnp(ik))
  enddo
  close(160)
endif

! find maximum number of local orbitals over all l-channels
allocate(ltmp(0:lolmax))
nlomaxl=0
do is=1,nspecies
  ltmp=0
  do ilo=1,nlorb(is)
    ltmp(lorbl(ilo,is))=ltmp(lorbl(ilo,is))+1
  enddo
  do l=0,lolmax
    if (l.le.lmaxvr) then
      nlomaxl=max(nlomaxl,ltmp(l))
    endif
  enddo
enddo
deallocate(ltmp)
! maximum number of radial functions for each l-channel
mtord=apwordmax+nlomaxl
if (iproc.eq.0) then
  write(150,*)
  write(150,'("Maximum number of l.o. over all l-channels : ",I4)')nlomaxl
  write(150,'("Maximum order of radial functions          : ",I4)')mtord
endif

allocate(uuj(0:lmaxvr,0:lmaxvr,0:lmaxexp,mtord,mtord,natmtot))

! find maximum number of Gaunt-like coefficients
ngumax=0
do l1=0,lmaxvr
  do m1=-l1,l1 
    lm1=idxlm(l1,m1)
    do l2=0,lmaxvr
      do m2=-l2,l2
        lm2=idxlm(l2,m2)
        do l3=0,lmaxexp
          do m3=-l3,l3
            lm3=idxlm(l3,m3)
	    if (abs(gaunt(l2,l1,l3,m2,m1,m3)).gt.1d-10) ngumax=ngumax+1
          enddo
        enddo
      enddo
    enddo
  enddo
enddo
ngumax=ngumax*mtord*mtord
if (iproc.eq.0) then
  write(150,*)
  write(150,'("Maximum number of Gaunt-like coefficients : ",I8)')ngumax
  call flushifc(150)
endif

allocate(ngu(natmtot))
allocate(gu(ngumax,natmtot))
allocate(igu(4,ngumax,natmtot))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
allocate(zrhofc0(max_num_nnp))
  
! different implementation for parallel and serial execution
#ifdef _MPI_

allocate(evecfvloc(nmatmax,nstfv,nspnfv,nkptlocnr(iproc)))
allocate(evecsvloc(nstsv,nstsv,nkptlocnr(iproc)))
allocate(acoeffloc(lmmaxvr,mtord,natmtot,nspinor,nstsv,nkptlocnr(iproc)))
allocate(acoeff2(lmmaxvr,mtord,natmtot,nspinor,nstsv))
allocate(status(MPI_STATUS_SIZE))
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
  call tcoeff1(ngknr(ikloc),mtord,apwalm,evecfvloc(1,1,1,ikloc), &
    evecsvloc(1,1,ikloc),acoeffloc(1,1,1,1,1,ikloc))
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

zrhofc=dcmplx(0.d0,0.d0)
do ig=1,ngvec_me
  if (iproc.eq.0) then
    write(150,*)
    write(150,'("G-vector ",I4," out of ",I4)')ig,ngvec_me
  endif
! calculate radial integrals <u|j|u'>
  call calc_uuj(uuj,lmaxexp,gq0(ig),mtord)
! compute Gaunt-like coefficients
  gu=dcmplx(0.d0,0.d0)
  igu=0
  ngu=0
  do ias=1,natmtot
    do io1=1,mtord
    do io2=1,mtord
      do l1=0,lmaxvr
      do m1=-l1,l1 
        lm1=idxlm(l1,m1)
        do l2=0,lmaxvr
        do m2=-l2,l2
          lm2=idxlm(l2,m2)
          do l3=0,lmaxexp
          do m3=-l3,l3
            lm3=idxlm(l3,m3)
	    t1=gaunt(l2,l1,l3,m2,m1,m3)*uuj(l1,l2,l3,io1,io2,ias)
            if (abs(t1).gt.1d-10) then
              ngu(ias)=ngu(ias)+1
              gu(ngu(ias),ias)=t1*ylmgq0(lm3,ig)*dconjg(zi**l3)*fourpi*dconjg(sfacgq0(ig,ias))
              igu(1,ngu(ias),ias)=lm1
              igu(2,ngu(ias),ias)=lm2
              igu(3,ngu(ias),ias)=io1
              igu(4,ngu(ias),ias)=io2
            endif
          enddo
          enddo
        enddo
        enddo
      enddo
      enddo
    enddo
    enddo
  enddo

  if (iproc.eq.0) then
    write(150,*)
    write(150,'("  Starting k-point loop")')
  endif
  
  do ikstep=1,nkptlocnr(0)
    if (iproc.eq.0) then
      write(150,'("    k-step ",I4," out of ",I4)')ikstep,nkptlocnr(0)
      call flushifc(150)
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
        call mpi_isend(acoeffloc(1,1,1,1,1,ik),                &
          lmmaxvr*mtord*natmtot*nspinor*nstsv,                 &
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
        call mpi_recv(acoeff2,lmmaxvr*mtord*natmtot*nspinor*nstsv,     &
          MPI_DOUBLE_COMPLEX,isend(ikstep,iproc,1),tag,MPI_COMM_WORLD,status,ierr)
        tag=tag+1
        call mpi_recv(ngknr2,1,MPI_INTEGER,isend(ikstep,iproc,1),tag,MPI_COMM_WORLD,status,ierr)
        tag=tag+1
        call mpi_recv(igkignr2,ngkmax,MPI_INTEGER,isend(ikstep,iproc,1),tag,MPI_COMM_WORLD,status,ierr)
      else
        ik=isend(ikstep,iproc,2)
        evecfv2(:,:,:)=evecfvloc(:,:,:,ik)
        evecsv2(:,:)=evecsvloc(:,:,ik)
        acoeff2(:,:,:,:,:)=acoeffloc(:,:,:,:,:,ik)
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
      
      zrhofc0=dcmplx(0.d0,0.d0)    
      call cpu_time(cpu0)
! calculate interstitial contribution for all combinations of n,n'
      call zrhoftistl_v3(ig,max_num_nnp,num_nnp(ikstep),nnp(1,1,ikstep),ngknr(ikstep),ngknr2, &
        igkignr(1,ikstep),igkignr2,ikq(ik,4),evecfvloc(1,1,1,ikstep),evecsvloc(1,1,ikstep), &
        evecfv2,evecsv2,zrhofc0)
      call cpu_time(cpu1)
      timeistl=cpu1-cpu0
    
      call cpu_time(cpu0)
      call zrhoftmt_v3(ngvec_me,max_num_nnp,num_nnp(ikstep),nnp(1,1,ikstep),mtord, &
        acoeffloc(1,1,1,1,1,ikstep),acoeff2,ngumax,ngu,gu,igu,zrhofc0)
      call cpu_time(cpu1)
      timemt=cpu1-cpu0
      
      zrhofc(ig,:,ikstep)=zrhofc0(:)

      if (iproc.eq.0) then
        write(150,'("  interstitial time (seconds) : ",F12.2)')timeistl
        write(150,'("    muffin-tin time (seconds) : ",F12.2)')timemt
        call flushifc(150)
      endif
    
    endif ! (ikstep.le.nkptlocnr(iproc))
  
  enddo !ikstep
enddo !ig
call pstop

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
deallocate(acoeffloc)
deallocate(status)
deallocate(acoeff2)
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
allocate(acoeff1(lmmaxvr,mtord,natmtot,nspinor,nstsv))
allocate(acoeff2(lmmaxvr,mtord,natmtot,nspinor,nstsv))

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
  call tcoeff1(ngknr(ik),mtord,apwalm,evecfv1,evecsv1,acoeff1)

! generate data at k'=k+q-K
  call getevecfv(vklnr(1,jk),vgklnr(1,1,jk),evecfv2)
  call getevecsv(vklnr(1,jk),evecsv2) 
  call match(ngknr(jk),gknr(1,jk),tpgknr(1,1,jk),sfacgknr(1,1,jk),apwalm)
  call tcoeff1(ngknr(jk),mtord,apwalm,evecfv2,evecsv2,acoeff2)
   
  zrhofc=dcmplx(0.d0,0.d0)
  call cpu_time(cpu0)
! calculate interstitial contribution for all combinations of n,n'
  call zrhoftistl(ngvec_me,max_num_nnp,num_nnp(ik),nnp(1,1,ik),ngknr(ik),ngknr(jk), &
    igkignr(1,ik),igkignr(1,jk),ikq(ik,4),evecfv1,evecsv1,evecfv2,evecsv2,zrhofc)
  call cpu_time(cpu1)
  timeistl=cpu1-cpu0

  call cpu_time(cpu0)  
! calculate MT contribution fot all combinations of n,n'
!  call zrhoftmt1(ngvec_me,max_num_nnp,num_nnp(ik),nnp(1,1,ik),mtord,ylmgq0, &
!    sfacgq0,acoeff1,acoeff2,uuj,gnt,lgnt,zrhofc)
  call zrhoftmt2(ngvec_me,max_num_nnp,num_nnp(ik),nnp(1,1,ik),mtord,ylmgq0, &
    sfacgq0,acoeff1,acoeff2,uuj,ngumax,ngu,gu,igu,zrhofc)
  call cpu_time(cpu1)
  timemt=cpu1-cpu0
  
  write(150,'("  interstitial time (seconds) : ",F12.2)')timeistl
  write(150,'("    muffin-tin time (seconds) : ",F12.2)')timemt
  
  write(160)zrhofc(1:ngvec_me,1:num_nnp(ik),1)
enddo !ik
close(160)

deallocate(evecfv1)
deallocate(evecsv1)
deallocate(evecfv2)
deallocate(evecsv2)
deallocate(zrhofc)
deallocate(acoeff1)
deallocate(acoeff2)

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


subroutine tcoeff1(ngp,mtord,apwalm,evecfv,evecsv,acoeff)
use modmain
implicit none
! arguments
integer, intent(in) :: ngp
integer, intent(in) :: mtord
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
complex(8), intent(out) :: acoeff(lmmaxvr,mtord,natmtot,nspinor,nstsv)
! local variables
integer j,l,m,ispn,istfv,is,ia,ias,lm,ig,i1,io,ilo
integer ordl(0:lmaxvr)
complex(8), allocatable :: acoeff_t(:,:,:,:)

allocate(acoeff_t(nstfv,mtord,lmmaxvr,natmtot))
acoeff_t=dcmplx(0.d0,0.d0)
! calculate first-variational coefficients
do istfv=1,nstfv
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      ordl=0
! apw coefficients
      do l=0,lmaxvr
        do io=1,apword(l,is)
          ordl(l)=ordl(l)+1
          do m=-l,l
            lm=idxlm(l,m)
            do ig=1,ngp
              acoeff_t(istfv,ordl(l),lm,ias)=acoeff_t(istfv,ordl(l),lm,ias) + &
                evecfv(ig,istfv)*apwalm(ig,io,lm,ias)
            enddo !ig
          enddo !m
        enddo !io
      enddo !l
! local orbital coefficients     
      do ilo=1,nlorb(is)
        l=lorbl(ilo,is)
        if (l.le.lmaxvr) then
          ordl(l)=ordl(l)+1
          do m=-l,l
            lm=idxlm(l,m)
            i1=ngp+idxlo(lm,ilo,ias)
            acoeff_t(istfv,ordl(l),lm,ias)=evecfv(i1,istfv)
          enddo !m
        endif
      enddo !ilo
    enddo !ia 
  enddo !is
enddo !istfv
! calculate second-variational coefficients
acoeff=dcmplx(0.d0,0.d0)
do j=1,nstsv
  do ispn=1,nspinor
    do ias=1,natmtot
      do io=1,mtord
        do lm=1,lmmaxvr
          do istfv=1,nstfv
	    acoeff(lm,io,ias,ispn,j)=acoeff(lm,io,ias,ispn,j) + &
	      evecsv(istfv+(ispn-1)*nstfv,j)*acoeff_t(istfv,io,lm,ias)
	  enddo !istfv
	enddo !lm
      enddo !io
    enddo !ias
  enddo !ispn
enddo !j
deallocate(acoeff_t)
return
end


subroutine calc_uuj(uuj,lmaxexp,gq0,mtord)
use modmain
implicit none
! arguments
integer, intent(in) :: lmaxexp
integer, intent(in) :: mtord
real(8), intent(in) :: gq0
real(8), intent(out) :: uuj(0:lmaxvr,0:lmaxvr,0:lmaxexp,mtord,mtord,natmtot)
! local variables
integer ia,is,ias,l,io,ilo,l1,l2,l3,io1,io2,ir
real(8), allocatable :: ufr(:,:,:)
real(8), allocatable :: jlgq0r(:,:,:)
integer ordl(0:lmaxvr)
real(8) fr(nrmtmax),gr(nrmtmax),cf(3,nrmtmax)
real(8) t1,jl(0:lmaxexp)

allocate(ufr(nrmtmax,0:lmaxvr,mtord))
allocate(jlgq0r(nrmtmax,0:lmaxexp,nspecies))

! generate Bessel functions j_l(|G+q'|x)
do is=1,nspecies
  do ir=1,nrmt(is)
    t1=gq0*spr(ir,is)
    call sbessel(lmaxexp,t1,jl)
    jlgq0r(ir,:,is)=jl(:)
  enddo
enddo

uuj=0.d0

do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! store all radial functions in one array
    ufr=0.d0
    ordl=0
! apw functions
    do l=0,lmaxvr
      do io=1,apword(l,is)
        ordl(l)=ordl(l)+1
        ufr(1:nrmt(is),l,ordl(l))=apwfr(1:nrmt(is),1,io,l,ias)
      enddo
    enddo
! lo functions
    do ilo=1,nlorb(is)
      l=lorbl(ilo,is)
      if (l.le.lmaxvr) then
        ordl(l)=ordl(l)+1
        ufr(1:nrmt(is),l,ordl(l))=lofr(1:nrmt(is),1,ilo,ias)
      endif
    enddo
    do l1=0,lmaxvr
      do l2=0,lmaxvr
        do l3=0,lmaxexp
          do io1=1,mtord
            do io2=1,mtord
              do ir=1,nrmt(is)
                fr(ir)=ufr(ir,l1,io1)*ufr(ir,l2,io2)*jlgq0r(ir,l3,is)*(spr(ir,is)**2)
              enddo
              call fderiv(-1,nrmt(is),spr(1,is),fr,gr,cf)
              uuj(l1,l2,l3,io1,io2,ias)=gr(nrmt(is))
            enddo !io2
          enddo !io1
        enddo !l3
      enddo !l2
    enddo !l1   
  enddo !ia
enddo !is

deallocate(ufr,jlgq0r)
return
end   

    
subroutine zrhoftmt_v3(ngvec_me,max_num_nnp,num_nnp,nnp,mtord, &
  acoeff1,acoeff2,ngumax,ngu,gu,igu,zrhofc)
use modmain
implicit none
! arguments
integer, intent(in) :: ngvec_me
integer, intent(in) :: max_num_nnp
integer, intent(in) :: num_nnp
integer, intent(in) :: nnp(max_num_nnp,3)
integer, intent(in) :: mtord
integer, intent(in) :: ngumax
integer, intent(in) :: ngu(natmtot)
integer, intent(in) :: igu(4,ngumax,natmtot)
complex(8), intent(in) :: gu(ngumax,natmtot)
complex(8), intent(in) :: acoeff1(lmmaxvr,mtord,natmtot,nspinor,nstsv)
complex(8), intent(in) :: acoeff2(lmmaxvr,mtord,natmtot,nspinor,nstsv)
complex(8), intent(inout) :: zrhofc(max_num_nnp)
! local variables
integer i,j,ist1,ist2,isp,ias,io1,io2,l1,m1,lm1,l2,m2,lm2,l3,m3,lm3,ispn 
complex(8) zt1,zt2,zt3,zt4
complex(8) a1(lmmaxvr,mtord),a2(lmmaxvr,mtord)

do i=1,num_nnp
  ist1=nnp(i,1)
  ist2=nnp(i,2)
  do ispn=1,nspinor
    do ias=1,natmtot
      a1=dconjg(acoeff1(:,:,ias,ispn,ist1))
      a2=acoeff2(:,:,ias,ispn,ist2)
      do j=1,ngu(ias)
        lm1=igu(1,j,ias)
        lm2=igu(2,j,ias)
        io1=igu(3,j,ias)
        io2=igu(4,j,ias)
        zrhofc(i)=zrhofc(i)+a1(lm1,io1)*a2(lm2,io2)*gu(j,ias)
      enddo
    enddo !ias 
  enddo !ispn
enddo !i

return
end     

subroutine zrhoftistl_v3(ig,max_num_nnp,num_nnp,nnp,ngknri,ngknrj,igkignri,igkignrj, &
  igkq,evecfvi,evecsvi,evecfvj,evecsvj,zrhofc)
use modmain
implicit none
! arguments
integer, intent(in) :: ig
integer, intent(in) :: max_num_nnp
integer, intent(in) :: num_nnp
integer, intent(in) :: nnp(max_num_nnp,3)
integer, intent(in) :: ngknri
integer, intent(in) :: ngknrj
integer, intent(in) :: igkq
integer, intent(in) :: igkignri(ngkmax)
integer, intent(in) :: igkignrj(ngkmax)
complex(8), intent(in) :: evecfvi(nmatmax,nstfv)
complex(8), intent(in) :: evecsvi(nstsv,nstsv)
complex(8), intent(in) :: evecfvj(nmatmax,nstfv)
complex(8), intent(in) :: evecsvj(nstsv,nstsv)
complex(8), intent(inout) :: zrhofc(max_num_nnp)

complex(8), allocatable :: mit(:,:)
complex(8), allocatable :: mit1(:,:)

integer is,ia,ias,igi,igj,ist1,ist2,i,i1,i2,ispn
integer iv4g(3)
real(8) v1(3),v2(3),tp4g(2),len4g
complex(8) sfac4g(natmtot)
complex(8) zt1

allocate(mit(ngknri,ngknrj))
allocate(mit1(nstfv,nstfv))
mit=dcmplx(0.d0,0.d0)
do igi=1,ngknri
  do igj=1,ngknrj
    ! G1-G2+G+K
    iv4g(:)=ivg(:,igkignri(igi))-ivg(:,igkignrj(igj))+ivg(:,ig)+ivg(:,igkq)
    if (sum(abs(iv4g)).eq.0) mit(igi,igj)=dcmplx(1.d0,0.d0)
    v2(:)=1.d0*iv4g(:)
    call r3mv(bvec,v2,v1)
    call sphcrd(v1,len4g,tp4g)
    call gensfacgp(1,v1,1,sfac4g)
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
	if (len4g.lt.1d-8) then
	  mit(igi,igj)=mit(igi,igj)-(fourpi/omega)*dconjg(sfac4g(ias))*(rmt(is)**3)/3.d0
	else
	  mit(igi,igj)=mit(igi,igj)-(fourpi/omega)*dconjg(sfac4g(ias)) * &
	    (-(rmt(is)/len4g**2)*cos(len4g*rmt(is))+(1/len4g**3)*sin(len4g*rmt(is)))
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
        zt1=zt1+dconjg(evecsvi(i1+(ispn-1)*nstfv,ist1))*mit1(i1,i2)
      enddo
      zrhofc(i)=zrhofc(i)+zt1*evecsvj(i2+(ispn-1)*nstfv,ist2)
    enddo
  enddo
enddo

deallocate(mit,mit1)
return
end        
  
