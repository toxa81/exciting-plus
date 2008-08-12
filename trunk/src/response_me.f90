subroutine response_me(ivq0m,ngvec_me)
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
integer, allocatable :: ikq(:,:)
! indexes for fft-transform of u_{nk}^{*}u_{n'k'}exp{-iKx}, where k'=k+q-K
integer, allocatable :: igfft1(:,:)
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
real(8), allocatable :: vgklnr(:,:,:)
real(8), allocatable :: vgkcnr(:,:)
real(8), allocatable :: gknr(:,:)
real(8), allocatable :: tpgknr(:,:,:)
real(8), allocatable :: occsvnr(:,:)
real(8), allocatable :: docc(:,:)
complex(8), allocatable :: sfacgknr(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfmt1(:,:,:,:,:)
complex(8), allocatable :: wfmt2(:,:,:,:,:)
complex(8), allocatable :: wfir1(:,:,:)
complex(8), allocatable :: wfir2(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)
complex(8), allocatable :: zrhofc(:,:,:)
complex(8), allocatable :: zrhofc1(:,:)
complex(8), allocatable :: evecfv0(:,:,:,:)
complex(8), allocatable :: evecsv0(:,:,:)
complex(8), allocatable :: evecfv1(:,:,:)
complex(8), allocatable :: evecsv1(:,:)
complex(8), allocatable :: evecfv2(:,:,:)
complex(8), allocatable :: evecsv2(:,:)

integer i,j,i1,ik,jk,ig,is,ir,ikstep,ist1,ist2,ispn
real(8) vkq0l(3),t1,jl(0:lmaxvr)
integer ivg1(3),ivg2(3)
complex(8) znorm

! for parallel execution
integer, allocatable :: isend(:,:,:)
integer tag,req,ierr
integer, allocatable :: status(:)
character*100 :: fname

! external functions
complex(8), external :: zfint
real(8), external :: r3taxi

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

! find number of G-vectors by given number of G-shells
!call getngvec(ngsh_me,ngvec_me)

! generate G+k vectors for entire BZ (this is required to compute 
!   wave-functions at each k-point)
allocate(vgklnr(3,ngkmax,nkptnr))
allocate(vgkcnr(3,ngkmax))
allocate(gknr(ngkmax,nkptnr))
allocate(tpgknr(2,ngkmax,nkptnr))
allocate(ngknr(nkptnr))
allocate(sfacgknr(ngkmax,natmtot,nkptnr))
allocate(igkignr(ngkmax,nkptnr))
do ik=1,nkptnr
  call gengpvec(vklnr(1,ik),vkcnr(1,ik),ngknr(ik),igkignr(1,ik), &
    vgklnr(1,1,ik),vgkcnr,gknr(1,ik),tpgknr(1,1,ik))
  call gensfacgp(ngknr(ik),vgkcnr,ngkmax,sfacgknr(1,1,ik))
enddo

! allocate memory for wafe-functions and complex charge density
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngrtot,nspinor,nstsv))
allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir2(ngrtot,nspinor,nstsv))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
allocate(zrhoir(ngrtot))
  
allocate(ikq(nkptnr,3))
allocate(vgq0c(3,ngvec))
allocate(gq0(ngvec))
allocate(tpgq0(2,ngvec))
allocate(sfacgq0(ngvec,natmtot))
allocate(ylmgq0(lmmaxvr,ngvec)) 
allocate(jlgq0r(nrcmtmax,0:lmaxvr,nspecies,ngvec))
allocate(occsvnr(nstsv,nkptnr))
allocate(igfft1(ngvec_me,nkptnr))
allocate(zrhofc1(ngvec_me,3))

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
! search for new fft indexes
  do ig=1,ngvec_me
    ivg2(:)=ivg(:,ig)+ivg1(:)
    igfft1(ig,ik)=igfft(ivgig(ivg2(1),ivg2(2),ivg2(3)))
  enddo
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
allocate(nnp(nkptnr,max_num_nnp,3))
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
          nnp(ik,i1,1)=ist1
          nnp(ik,i1,2)=ist2
	  if (spin_me.eq.3) then
	    nnp(ik,i1,3)=ispn
	  else
            nnp(ik,i1,3)=1
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
    write(160)nnp(ik,1:num_nnp(ik),1:3)
    write(160)docc(ik,1:num_nnp(ik))
  enddo
  close(160)
endif
  
! different implementation for parallel and serial execution
#ifdef _MPI_

allocate(nkptloc(0:nproc-1))
allocate(ikptloc(0:nproc-1,2))
call splitk(nkptnr,nproc,nkptloc,ikptloc)
    
! root proc reads all eigen-vectors to memory
if (iproc.eq.0) then
  allocate(evecfv0(nmatmax,nstfv,nspnfv,nkpt))
  allocate(evecsv0(nstsv,nstsv,nkpt))
  do ik=1,nkpt
    call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv0(1,1,1,ik))
    call getevecsv(vkl(1,ik),evecsv0(1,1,ik))
  enddo 
endif
    
! find indexes of k-points to send and receive
allocate(isend(nkptloc(0),0:nproc-1,2))
do ikstep=1,nkptloc(0)
  do i=0,nproc-1
    ik=ikptloc(i,1)+ikstep-1
    if (ikstep.le.nkptloc(i)) then
! for the ikstep the i-th proc requires eigen-vectors in this 
!   two irreducible k-points from root proc
      isend(ikstep,i,1)=ikq(ik,2)
      isend(ikstep,i,2)=ikq(ik,3)
    endif
  enddo
enddo
    
allocate(status(MPI_STATUS_SIZE))
allocate(evecfv1(nmatmax,nstfv,nspnfv))
allocate(evecsv1(nstsv,nstsv))
allocate(evecfv2(nmatmax,nstfv,nspnfv))
allocate(evecsv2(nstsv,nstsv))
allocate(zrhofc(ngvec_me,max_num_nnp,nkptloc(iproc)))

if (iproc.eq.0) then
  write(150,*)
  write(150,'("Starting k-point loop")')
endif
do ikstep=1,nkptloc(0)
  if (iproc.eq.0) then
    write(150,'("k-step ",I4," out of ",I4)')ikstep,nkptloc(0)
    call flushifc(150)
! root proc should send eigen-vectors
    do i=1,nproc-1
      if (ikstep.le.nkptloc(i)) then
! send eigen-vectors at two k-points to proc i 
        tag=(ikstep*nproc+i)*10
	ik=isend(ikstep,i,1)
	call mpi_isend(evecfv0(1,1,1,ik),nmatmax*nstfv*nspnfv,MPI_DOUBLE_COMPLEX,i,tag,MPI_COMM_WORLD,req,ierr)
	tag=tag+1
	call mpi_isend(evecsv0(1,1,ik),nstsv*nstsv,MPI_DOUBLE_COMPLEX,i,tag,MPI_COMM_WORLD,req,ierr)
	tag=tag+1
	ik=isend(ikstep,i,2)
	call mpi_isend(evecfv0(1,1,1,ik),nmatmax*nstfv*nspnfv,MPI_DOUBLE_COMPLEX,i,tag,MPI_COMM_WORLD,req,ierr)
	tag=tag+1
	call mpi_isend(evecsv0(1,1,ik),nstsv*nstsv,MPI_DOUBLE_COMPLEX,i,tag,MPI_COMM_WORLD,req,ierr)
      endif
    enddo !i
! copy local arrays
    ik=isend(ikstep,0,1)
    evecfv1(:,:,:)=evecfv0(:,:,:,ik)
    evecsv1(:,:)=evecsv0(:,:,ik)
    ik=isend(ikstep,0,2)
    evecfv2(:,:,:)=evecfv0(:,:,:,ik)
    evecsv2(:,:)=evecsv0(:,:,ik)
  else
    if (ikstep.le.nkptloc(iproc)) then
      tag=(ikstep*nproc+iproc)*10
      call mpi_recv(evecfv1,nmatmax*nstfv*nspnfv,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,status,ierr)
      tag=tag+1
      call mpi_recv(evecsv1,nstsv*nstsv,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,status,ierr)
      tag=tag+1
      call mpi_recv(evecfv2,nmatmax*nstfv*nspnfv,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,status,ierr)
      tag=tag+1
      call mpi_recv(evecsv2,nstsv*nstsv,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,status,ierr)
    endif
  endif
  call mpi_barrier(MPI_COMM_WORLD,ierr)
      
  if (ikstep.le.nkptloc(iproc)) then
    ik=ikptloc(iproc,1)+ikstep-1
    jk=ikq(ik,1)

! generate wave-functions at k
    call getevecfvp(vklnr(1,ik),vgklnr(1,1,ik),evecfv1,ikq(ik,2))
    call getevecsvp(vklnr(1,ik),evecsv1) 
    call match(ngknr(ik),gknr(1,ik),tpgknr(1,1,ik),sfacgknr(1,1,ik),apwalm)
    call genwfsv(.false.,ngknr(ik),igkignr(1,ik),evalsv(1,1),apwalm,evecfv1, &
      evecsv1,wfmt1,wfir1)

! test normalization    
    do i=1,nstsv
      call vnlrho(.true.,wfmt1(1,1,1,1,i),wfmt1(1,1,1,1,i),wfir1(1,1,i), &
        wfir1(1,1,i),zrhomt,zrhoir)
      znorm=zfint(zrhomt,zrhoir)
      if (abs(znorm-1.d0).gt.0.1d0) then
        write(*,*)
        write(*,'("Error(response_me): bad norm ",G18.10," of wave-function ",&
          & I4," at k-point ",I4)')abs(znorm),i,ik
	write(*,'("  try to switch off the symmetry")')
	write(*,*)
	call pstop
      endif
    enddo

! generate wave-functions at k'=k+q-K
    call getevecfvp(vklnr(1,jk),vgklnr(1,1,jk),evecfv2,ikq(ik,3))
    call getevecsvp(vklnr(1,jk),evecsv2) 
    call match(ngknr(jk),gknr(1,jk),tpgknr(1,1,jk),sfacgknr(1,1,jk),apwalm)
    call genwfsv(.false.,ngknr(jk),igkignr(1,jk),evalsv(1,1),apwalm,evecfv2, &
      evecsv2,wfmt2,wfir2)

    do i=1,num_nnp(ik)
      ist1=nnp(ik,i,1)
      ist2=nnp(ik,i,2)
      call vnlrho(.true.,wfmt1(1,1,1,1,ist1),wfmt2(1,1,1,1,ist2),wfir1(1,1,ist1), &
        wfir2(1,1,ist2),zrhomt,zrhoir)
      call zrhoft(zrhomt,zrhoir,jlgq0r,ylmgq0,sfacgq0,ngvec_me,igfft1(1,ik),zrhofc1)
      zrhofc(:,i,ikstep)=zrhofc1(:,3)
    enddo

  endif ! (ikstep.le.nkptloc(iproc))
enddo !ikstep

call mpi_barrier(MPI_COMM_WORLD,ierr)

do i=0,nproc-1
  if (i.eq.iproc) then
    open(160,file=trim(fname),form='unformatted',status='old',position='append')
    do ikstep=1,nkptloc(iproc)
      ik=ikptloc(iproc,1)+ikstep-1
      write(160)ik,ikq(ik,1)
      write(160)zrhofc(1:ngvec_me,1:num_nnp(ik),ikstep)
    enddo !ikstep
    close(160)
  endif !i.eq.iproc
  call mpi_barrier(MPI_COMM_WORLD,ierr)
enddo 
  
if (iproc.eq.0) then
  deallocate(evecfv0)
  deallocate(evecsv0)
endif 
deallocate(nkptloc)
deallocate(ikptloc)
deallocate(status)
deallocate(isend)
deallocate(evecfv1)
deallocate(evecsv1)
deallocate(evecfv2)
deallocate(evecsv2)
deallocate(zrhofc)

#else

allocate(evecfv1(nmatmax,nstfv,nspnfv))
allocate(evecsv1(nstsv,nstsv))
allocate(zrhofc(ngvec_me,max_num_nnp,1))

open(160,file=trim(fname),form='unformatted',status='old',position='append')

write(150,*)
do ik=1,nkptnr
  write(150,'("k-point ",I4," out of ",I4)')ik,nkptnr
  
  jk=ikq(ik,1)
  
  write(160)ik,jk
  
! generate wave-functions at k
  call getevecfv(vklnr(1,ik),vgklnr(1,1,ik),evecfv1)
  call getevecsv(vklnr(1,ik),evecsv1) 
  call match(ngknr(ik),gknr(1,ik),tpgknr(1,1,ik),sfacgknr(1,1,ik),apwalm)
  call genwfsv(.false.,ngknr(ik),igkignr(1,ik),evalsv(1,1),apwalm,evecfv1, &
    evecsv1,wfmt1,wfir1)

! test normalization    
  do i=1,nstsv
    call vnlrho(.true.,wfmt1(1,1,1,1,i),wfmt1(1,1,1,1,i),wfir1(1,1,i), &
      wfir1(1,1,i),zrhomt,zrhoir)
    znorm=zfint(zrhomt,zrhoir)
    if (abs(znorm-1.d0).gt.0.01d0) then
      write(150,'("Warning: bad norm ",G18.10," of wave-function ",&
        & I4," at k-point ",I4)')abs(znorm),i,ik
    endif
  enddo

! generate wave-functions at k'=k+q-K
  call getevecfv(vklnr(1,jk),vgklnr(1,1,jk),evecfv1)
  call getevecsv(vklnr(1,jk),evecsv1) 
  call match(ngknr(jk),gknr(1,jk),tpgknr(1,1,jk),sfacgknr(1,1,jk),apwalm)
  call genwfsv(.false.,ngknr(jk),igkignr(1,jk),evalsv(1,1),apwalm,evecfv1, &
    evecsv1,wfmt2,wfir2)
  
  do i=1,num_nnp(ik)
    ist1=nnp(ik,i,1)
    ist2=nnp(ik,i,2)
    call vnlrho(.true.,wfmt1(1,1,1,1,ist1),wfmt2(1,1,1,1,ist2),wfir1(1,1,ist1), &
      wfir2(1,1,ist2),zrhomt,zrhoir)
    call zrhoft(zrhomt,zrhoir,jlgq0r,ylmgq0,sfacgq0,ngvec_me,igfft1(1,ik),zrhofc1)
    zrhofc(:,i,1)=zrhofc1(:,3)
  enddo
  
  write(160)zrhofc(1:ngvec_me,1:num_nnp(ik),1)
enddo !ik
close(160)
deallocate(evecfv1)
deallocate(evecsv1)
deallocate(zrhofc)

#endif

deallocate(vgklnr)
deallocate(vgkcnr)
deallocate(gknr)
deallocate(tpgknr)
deallocate(ngknr)
deallocate(sfacgknr)
deallocate(igkignr)
deallocate(apwalm)
deallocate(wfmt1)
deallocate(wfir1)
deallocate(wfmt2)
deallocate(wfir2)
deallocate(zrhomt)
deallocate(zrhoir)
deallocate(ikq)
deallocate(vgq0c)
deallocate(gq0)
deallocate(tpgq0)
deallocate(sfacgq0)
deallocate(ylmgq0) 
deallocate(jlgq0r)
deallocate(occsvnr)
deallocate(igfft1)
deallocate(zrhofc1)
deallocate(num_nnp)
deallocate(nnp)
deallocate(docc)

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


subroutine getgshells(ngsh,igishell,ishellng)
use modmain
implicit none
! arguments
integer, intent(out) :: ngsh
integer, intent(out) :: igishell(ngvec)
integer, intent(out) :: ishellng(ngvec,2)

integer ish,ig

ish=1
ig=0
do while (ig.lt.ngvec)
  ig=ig+1
  igishell(ig)=ish
  if (abs(gc(ig+1)-gc(ig)).gt.epslat) then
    ish=ish+1
  endif
enddo 

ngsh=ish-1

do ish=1,ngsh
  ishellng(ish,1)=0
  do ig=1,ngvec
    if (igishell(ig).eq.ish) ishellng(ish,1)=ishellng(ish,1)+1
  enddo
  ishellng(ish,2)=sum(ishellng(1:ish,1))
enddo
 
return
end


subroutine getngvec(ngsh,ngv)
use modmain
implicit none
! arguments
integer ,intent(in)  :: ngsh
integer ,intent(out) :: ngv
! local variables
integer ,allocatable :: gshell(:)
integer              :: i

allocate(gshell(ngvec))

i=1
ngv=0
do while (i.le.ngsh)
  ngv=ngv+1
  gshell(ngv)=i
  if (abs(gc(ngv+1)-gc(ngv)).gt.epslat) then
    i=i+1
  endif
enddo 

if (iproc.eq.0) then
  write(150,*)
  write(150,'("Number of G-shells  : ",I4)')ngsh
  write(150,'("Number of G-vectors : ",I4)')ngv
  write(150,*)
  write(150,'("  G-vec.       lat.coord.      length(a.u.)   shell ")')
  write(150,'(" ---------------------------------------------------")')
  do i=1,ngv
    write(150,'(2X,I4,4X,3I5,4X,F12.6,5x,I4)')i,ivg(:,i),gc(i),gshell(i)
  enddo
endif

deallocate(gshell)

return
end

subroutine pstop
#ifdef _MPI_
use mpi
#endif
implicit none
integer ierr

#ifdef _MPI_
call mpi_barrier(MPI_COMM_WORLD,ierr)
call mpi_abort(MPI_COMM_WORLD,-1,ierr)
call mpi_finalize(ierr)
#else
stop
#endif

return
end

subroutine splitk(nkpt,nproc,nkptloc,ikptloc)
implicit none
! arguments
integer, intent(in) :: nkpt
integer, intent(in) :: nproc
integer, intent(out) :: nkptloc(0:nproc-1)
integer, intent(out) :: ikptloc(0:nproc-1,2) 
! local variables
integer i,n1,n2

! minimum number of k-points for each proc.
n1=nkpt/nproc
! remaining number of k-points which will be distributed among first n2 procs
n2=nkpt-n1*nproc
! each proc gets n1 k-points
nkptloc(:)=n1
! additionally, first n2 procs get extra point
do i=0,n2-1
  nkptloc(i)=nkptloc(i)+1
enddo 
! build index of first and last k-point for each proc
ikptloc(0,1)=1
ikptloc(0,2)=nkptloc(0)
do i=1,nproc-1
  ikptloc(i,1)=ikptloc(i-1,2)+1
  ikptloc(i,2)=ikptloc(i,1)+nkptloc(i)-1
enddo
    
!if (iproc.eq.0) then
!  write(150,*)
!  write(150,'(" iproc  first k   last k   nkpt ")')
!  write(150,'(" ------------------------------ ")')
!  do i=0,nproc-1
!    write(150,'(1X,I4,4X,I4,5X,I4,5X,I4)')i,ikptloc(i,1),ikptloc(i,2),nkptloc(i)
!  enddo
!endif

return
end
