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
!              2: index of K-vector which brings k+q to first BZ
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
complex(8), allocatable :: evecfv1(:,:,:)
complex(8), allocatable :: evecsv1(:,:)
complex(8), allocatable :: evecfv2(:,:,:)
complex(8), allocatable :: evecsv2(:,:)
complex(8), allocatable :: evecloc(:,:,:,:)
complex(8), allocatable :: evec2(:,:,:)
complex(8), allocatable :: acoeffloc(:,:,:,:,:,:)
complex(8), allocatable :: acoeff1(:,:,:,:,:)
complex(8), allocatable :: acoeff2(:,:,:,:,:)

integer i,j,i1,ik,jk,ig,is,ir,ikstep,ist1,ist2,ispn,ikloc,ilo,l,ia,ias,io1,io2,istfv
integer ngknr2
real(8) vkq0l(3),t1,jl(0:lmaxvr)
integer ivg1(3),ivg2(3)
real(8) cpu0,cpu1,timeistl,timemt
integer nlomaxl
integer mtord
integer, allocatable :: ltmp(:)
real(8), allocatable :: uuj(:,:,:,:,:,:,:)
complex(4), allocatable :: gu(:,:,:)
integer, allocatable :: igu(:,:,:,:)
integer, allocatable :: ngu(:,:)
integer ngumax

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
  write(150,*)
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

allocate(ikq(nkptnr,2))
allocate(vgq0c(3,ngvec_me))
allocate(gq0(ngvec_me))
allocate(tpgq0(2,ngvec_me))
allocate(sfacgq0(ngvec_me,natmtot))
allocate(ylmgq0(lmmaxexp,ngvec_me)) 
allocate(occsvnr(nstsv,nkptnr))
allocate(igfft1(ngvec_me,nkptnr))

! q-vector in lattice coordinates
do i=1,3
  vq0l(i)=1.d0*ivq0m(i)/ngridk(i)+1d-12
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
! search for new fft indexes
  do ig=1,ngvec_me
    ivg2(:)=ivg(:,ig)+ivg1(:)
    igfft1(ig,ik)=igfft(ivgig(ivg2(1),ivg2(2),ivg2(3)))
  enddo
  ikq(ik,2)=ivgig(ivg1(1),ivg1(2),ivg1(3))
!  if (iproc.eq.0) then
!    write(150,'(I4,2X,3F6.2,2X,3F6.2,2X,3I4,2X,3F6.2,2X,I4)') &
!    ik,vklnr(:,ik),vkq0l+ivg1,ivg1,vkq0l,ikq(ik,1)
!  endif
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
call getnnp(.true.,ikq,occsvnr,max_num_nnp,num_nnp,nnp,docc,&
  nkptlocnr(iproc),ikptlocnr(iproc,1))
#ifdef _MPI_
call mpi_allreduce(max_num_nnp,i,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
max_num_nnp=i
#endif  
allocate(num_nnp(nkptlocnr(iproc)))
allocate(nnp(max_num_nnp,3,nkptlocnr(iproc)))
allocate(docc(max_num_nnp,nkptlocnr(iproc)))
call getnnp(.false.,ikq,occsvnr,max_num_nnp,num_nnp,nnp,docc,&
  nkptlocnr(iproc),ikptlocnr(iproc,1))
if (iproc.eq.0) then
  write(150,*)
  write(150,'("Maximum number of interband transitions: ",I5)')max_num_nnp
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
  write(160)nkptnr,ngsh_me,ngvec_me,max_num_nnp,igq0
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

! find maximum number of local orbitals over all l-channels
call getmtord(lmaxvr,mtord)
!allocate(ltmp(0:lolmax))
!nlomaxl=0
!do is=1,nspecies
!  ltmp=0
!  do ilo=1,nlorb(is)
!    ltmp(lorbl(ilo,is))=ltmp(lorbl(ilo,is))+1
!  enddo
!  do l=0,lolmax
!    if (l.le.lmaxvr) then
!      nlomaxl=max(nlomaxl,ltmp(l))
!    endif
!  enddo
!enddo
!deallocate(ltmp)
!mtord=apwordmax+nlomaxl
if (iproc.eq.0) then
  write(150,*)
  write(150,'("Calculating radial integrals")')
  write(150,'("  maximum order of radial functions : ",I4)')mtord
endif
allocate(uuj(0:lmaxvr,0:lmaxvr,0:lmaxexp,mtord,mtord,natmtot,ngvec_me))
call calc_uuj(uuj,lmaxexp,gq0,mtord,ngvec_me)
if (iproc.eq.0) then
  write(150,'("Done.")')
  call flushifc(150)
endif

call getgu(.true.,ngvec_me,lmaxexp,mtord,uuj,ylmgq0,sfacgq0,ngumax,ngu,gu,igu)
allocate(ngu(natmtot,ngvec_me))
allocate(gu(ngumax,natmtot,ngvec_me))
allocate(igu(4,ngumax,natmtot,ngvec_me))
call getgu(.false.,ngvec_me,lmaxexp,mtord,uuj,ylmgq0,sfacgq0,ngumax,ngu,gu,igu)
if (iproc.eq.0) then
  write(150,*)
  write(150,'("Maximum number of Gaunt-like coefficients : ",I8)')ngumax
  call flushifc(150)
endif
deallocate(uuj)

allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
  
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
allocate(evecloc(nmatmax,nstsv,nspinor,nkptlocnr(iproc)))
allocate(evec2(nmatmax,nstsv,nspinor))


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
  call getacoeff(lmaxvr,lmmaxvr,ngknr(ikloc),mtord,apwalm,evecfvloc(1,1,1,ikloc), &
    evecsvloc(1,1,ikloc),acoeffloc(1,1,1,1,1,ikloc))
! hack to switch off l-channel
!  do is=1,nspecies
!    if (trim(adjustl(spsymb(is))).eq.'O') then
!      do ia=1,natoms(is)
!        ias=idxas(ia,is)
!	  acoeffloc(2:4,:,ias,:,29:55,ikloc)=dcmplx(0.d0,0.d0)
!      enddo
!    endif
!  enddo
  ias=7
  acoeffloc(2:4,:,ias,:,29:55,ikloc)=dcmplx(0.d0,0.d0)
  ias=8
  acoeffloc(2:4,:,ias,:,29:55,ikloc)=dcmplx(0.d0,0.d0)
enddo

if (iproc.eq.0) then
  write(150,'("Done.")')
  call flushifc(150)
endif

! transform eigen-vectors
if (iproc.eq.0) then
  write(150,*)
  write(150,'("Transforming eigen-vectors")')
  call flushifc(150)
endif
evecloc=dcmplx(0.d0,0.d0)
do ikloc=1,nkptlocnr(iproc)
  do ispn=1,nspinor                                                                                                                                
    do j=1,nstsv                                                                                                                                   
      do ig=1,ngknr(ikloc)                                                                                                                              
        do istfv=1,nstfv                                                                                                                           
          evecloc(ig,j,ispn,ikloc)=evecloc(ig,j,ispn,ikloc) + &
	    evecsvloc(istfv+(ispn-1)*nstfv,j,ikloc)*evecfvloc(ig,istfv,1,ikloc)                                               
        enddo
      enddo                                                                                                                                      
    enddo                                                                                                                                        
  enddo                                                                                                                                          
enddo
if (iproc.eq.0) then
  write(150,'("Done.")')
  call flushifc(150)
endif
deallocate(evecfvloc)
deallocate(evecsvloc)

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
zrhofc=dcmplx(0.d0,0.d0)
do ikstep=1,nkptlocnr(0)
  if (iproc.eq.0) then
    write(150,'("k-step ",I4," out of ",I4)')ikstep,nkptlocnr(0)
    call flushifc(150)
  endif
! find the i-th proc to which the current iproc should send data
  do i=0,nproc-1
    if (isend(ikstep,i,1).eq.iproc.and.iproc.ne.i) then
      tag=(ikstep*nproc+i)*10
      ik=isend(ikstep,i,2)
      call mpi_isend(evecloc(1,1,1,ik),                      &
        nmatmax*nstsv*nspinor,                               &
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
      call mpi_recv(evec2,nmatmax*nstsv*nspinor,MPI_DOUBLE_COMPLEX,  &
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
      evec2(:,:,:)=evecloc(:,:,:,ik)
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

! calculate interstitial contribution for all combinations of n,n'
    call cpu_time(cpu0)
    call zrhoftistl(ngvec_me,max_num_nnp,num_nnp(ikstep),nnp(1,1,ikstep),ngknr(ikstep),ngknr2, &
      igkignr(1,ikstep),igkignr2,ikq(ik,2),evecloc(1,1,1,ikstep),evec2,zrhofc(1,1,ikstep))
    call cpu_time(cpu1)
    timeistl=cpu1-cpu0

! calculate muffin-tin contribution for all combinations of n,n'    
    call cpu_time(cpu0)
    call zrhoftmt(ngvec_me,max_num_nnp,num_nnp(ikstep),nnp(1,1,ikstep),mtord, &
      acoeffloc(1,1,1,1,1,ikstep),acoeff2,ngumax,ngu,gu,igu,zrhofc(1,1,ikstep))
    call cpu_time(cpu1)
    timemt=cpu1-cpu0
  
    if (iproc.eq.0) then
      write(150,'("  interstitial time (seconds) : ",F12.2)')timeistl
      write(150,'("    muffin-tin time (seconds) : ",F12.2)')timemt
      call flushifc(150)
    endif
    
  endif ! (ikstep.le.nkptlocnr(iproc))
  
enddo !ikstep

call mpi_barrier(MPI_COMM_WORLD,ierr)

do i=0,nproc-1
  if (i.eq.iproc) then
    open(160,file=trim(fname),form='unformatted',status='old',position='append')
    do ikstep=1,nkptlocnr(iproc)
      ik=ikptlocnr(iproc,1)+ikstep-1
      write(160)ik,ikq(ik,1)
      write(160)num_nnp(ikstep)
      write(160)nnp(1:num_nnp(ikstep),1:3,ikstep)
      write(160)docc(1:num_nnp(ikstep),ikstep)
      write(160)zrhofc(1:ngvec_me,1:num_nnp(ikstep),ikstep)
    enddo !ikstep
    close(160)
  endif !i.eq.iproc
  call mpi_barrier(MPI_COMM_WORLD,ierr)
enddo 

deallocate(acoeffloc)
deallocate(acoeff2)
deallocate(status)
deallocate(evecfv2)
deallocate(evecsv2)
deallocate(igkignr2)
deallocate(zrhofc)
deallocate(evecloc)
deallocate(evec2)
  
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
  call getacoeff(lmaxvr,lmmaxvr,ngknr(ik),mtord,apwalm,evecfv1,evecsv1,acoeff1)

! generate data at k'=k+q-K
  call getevecfv(vklnr(1,jk),vgklnr(1,1,jk),evecfv2)
  call getevecsv(vklnr(1,jk),evecsv2) 
  call match(ngknr(jk),gknr(1,jk),tpgknr(1,1,jk),sfacgknr(1,1,jk),apwalm)
  call getacoeff(lmaxvr,lmmaxvr,ngknr(jk),mtord,apwalm,evecfv2,evecsv2,acoeff2)
   
  zrhofc=dcmplx(0.d0,0.d0)
  call cpu_time(cpu0)
! calculate interstitial contribution for all combinations of n,n'
  call zrhoftistl(ngvec_me,max_num_nnp,num_nnp(ik),nnp(1,1,ik),ngknr(ik),ngknr(jk), &
    igkignr(1,ik),igkignr(1,jk),ikq(ik,4),evecfv1,evecsv1,evecfv2,evecsv2,zrhofc,igfft1(1,ik))
  call cpu_time(cpu1)
  timeistl=cpu1-cpu0

  call cpu_time(cpu0)  
! calculate MT contribution fot all combinations of n,n'
!  call zrhoftmt1(ngvec_me,max_num_nnp,num_nnp(ik),nnp(1,1,ik),mtord,ylmgq0, &
!    sfacgq0,acoeff1,acoeff2,uuj,gnt,lgnt,zrhofc)
  call zrhoftmt(ngvec_me,max_num_nnp,num_nnp(ik),nnp(1,1,ik),mtord,ylmgq0, &
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


subroutine getacoeff(lmax,lmmax,ngp,mtord,apwalm,evecfv,evecsv,acoeff)
use modmain
implicit none
! arguments
integer, intent(in) :: lmax
integer, intent(in) :: lmmax
integer, intent(in) :: ngp
integer, intent(in) :: mtord
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
complex(8), intent(out) :: acoeff(lmmax,mtord,natmtot,nspinor,nstsv)
! local variables
integer j,l,m,ispn,istfv,is,ia,ias,lm,ig,i1,io,ilo
integer ordl(0:lmaxvr)
complex(8), allocatable :: acoeff_t(:,:,:,:)

allocate(acoeff_t(nstfv,mtord,lmmax,natmtot))
acoeff_t=dcmplx(0.d0,0.d0)
! calculate first-variational coefficients
do istfv=1,nstfv
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      ordl=0
! apw coefficients
      do l=0,lmax
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
        if (l.le.lmax) then
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
        do lm=1,lmmax
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

subroutine calc_uuj(uuj,lmaxexp,gq0,mtord,ngvec_me)
use modmain
implicit none
! arguments
integer, intent(in) :: lmaxexp
integer, intent(in) :: mtord
integer, intent(in) :: ngvec_me
real(8), intent(in) :: gq0(ngvec_me)
real(8), intent(out) :: uuj(0:lmaxvr,0:lmaxvr,0:lmaxexp,mtord,mtord,natmtot,ngvec_me)
! local variables
integer ia,is,ias,l,io,ilo,ig,l1,l2,l3,io1,io2,ir
real(8), allocatable :: ufr(:,:,:)
real(8), allocatable :: jlgq0r(:,:,:,:)
integer ordl(0:lmaxvr)
real(8) fr(nrmtmax),gr(nrmtmax),cf(3,nrmtmax)
real(8) t1,jl(0:lmaxexp)

allocate(ufr(nrmtmax,0:lmaxvr,mtord))
allocate(jlgq0r(nrmtmax,0:lmaxexp,nspecies,ngvec_me))

! generate Bessel functions j_l(|G+q'|x)
do ig=1,ngvec_me
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
    do ig=1,ngvec_me
      do l1=0,lmaxvr
        do l2=0,lmaxvr
          do l3=0,lmaxexp
            do io1=1,mtord
              do io2=1,mtord
                do ir=1,nrmt(is)
                  fr(ir)=ufr(ir,l1,io1)*ufr(ir,l2,io2)*jlgq0r(ir,l3,is,ig)*(spr(ir,is)**2)
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

deallocate(ufr,jlgq0r)
return
end   

subroutine zrhoftistl(ngvec_me,max_num_nnp,num_nnp,nnp,ngknri,ngknrj,igkignri,igkignrj, &
  ig1,eveci,evecj,zrhofc)
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
complex(8), intent(in) :: eveci(nmatmax,nstsv,nspinor)
complex(8), intent(in) :: evecj(nmatmax,nstsv,nspinor)
complex(8), intent(inout) :: zrhofc(ngvec_me,max_num_nnp)

complex(8), allocatable :: mit(:,:)
complex(8), allocatable :: a(:,:,:) 

integer is,ia,ias,ig,igi,igj,ist1,ist2,i,i1,i2,ispn,istfv
integer iv3g(3)
real(8) v1(3),v2(3),tp3g(2),len3g
complex(8) sfac3g(natmtot)
complex(8) zt1

allocate(a(ngknrj,nstsv,nspinor))
allocate(mit(ngknri,ngknrj))

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
  a=dcmplx(0.d0,0.d0)
  do ispn=1,nspinor
    do i=1,nstsv
      do igj=1,ngknrj
        do igi=1,ngknri
          a(igj,i,ispn)=a(igj,i,ispn)+dconjg(eveci(igi,i,ispn))*mit(igi,igj)
	enddo
      enddo
    enddo
  enddo
  do i=1,num_nnp
    ist1=nnp(i,1)
    ist2=nnp(i,2)
    do ispn=1,nspinor
      do igj=1,ngknrj
        zrhofc(ig,i)=zrhofc(ig,i)+evecj(igj,ist2,ispn)*a(igj,ist1,ispn)
      enddo
    enddo
  enddo
enddo
deallocate(mit,a)
return
end
        
subroutine zrhoftmt(ngvec_me,max_num_nnp,num_nnp,nnp,mtord, &
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
integer, intent(in) :: ngu(natmtot,ngvec_me)
integer, intent(in) :: igu(4,ngumax,natmtot,ngvec_me)
complex(4), intent(in) :: gu(ngumax,natmtot,ngvec_me)
complex(8), intent(in) :: acoeff1(lmmaxvr,mtord,natmtot,nspinor,nstsv)
complex(8), intent(in) :: acoeff2(lmmaxvr,mtord,natmtot,nspinor,nstsv)
complex(8), intent(inout) :: zrhofc(ngvec_me,max_num_nnp)
! local variables
integer ig,i,j,ist1,ist2,isp,ias,io1,io2,l1,m1,lm1,l2,m2,lm2,l3,m3,lm3,ispn 
complex(8) zt1,zt2,zt3,zt4
complex(8) a1(lmmaxvr,mtord),a2(lmmaxvr,mtord)

do ig=1,ngvec_me
  do i=1,num_nnp
    ist1=nnp(i,1)
    ist2=nnp(i,2)
    do ispn=1,nspinor
      do ias=1,natmtot
        a1=dconjg(acoeff1(:,:,ias,ispn,ist1))
        a2=acoeff2(:,:,ias,ispn,ist2)
        do j=1,ngu(ias,ig)
          lm1=igu(1,j,ias,ig)
          lm2=igu(2,j,ias,ig)
          io1=igu(3,j,ias,ig)
          io2=igu(4,j,ias,ig)
          zrhofc(ig,i)=zrhofc(ig,i)+a1(lm1,io1)*a2(lm2,io2)*gu(j,ias,ig)
        enddo
      enddo !ias 
    enddo !ispn
  enddo !i
enddo !ig    

return
end


subroutine getnnp(req,ikq,occsvnr,max_num_nnp,num_nnp,nnp,docc,nkptloc1,ikptloc1)
use modmain
implicit none
! arguments
logical, intent(in) :: req
integer, intent(in) :: ikq(nkptnr,2)
real(8), intent(in) :: occsvnr(nstsv,nkptnr)
integer, intent(in) :: nkptloc1
integer, intent(in) :: ikptloc1
integer, intent(inout) :: max_num_nnp
integer, intent(out) :: num_nnp(nkptnr)
integer, intent(out) :: nnp(max_num_nnp,3,nkptnr)
real(8), intent(out) :: docc(max_num_nnp,nkptnr)

integer i,ik,jk,istfv1,istfv2,ispn,ist1,ist2,ikloc

if (req) max_num_nnp=0
do ikloc=1,nkptloc1
  ik=ikptloc1+ikloc-1
  jk=ikq(ik,1)
  i=0
  do ispn=1,nspinor
    do istfv1=1,nstfv
      ist1=istfv1+(ispn-1)*nstfv
      do istfv2=1,nstfv
        ist2=istfv2+(ispn-1)*nstfv
        if ((ispn.eq.spin_me.or.spin_me.eq.3) .and. &
	    abs(occsvnr(ist1,ik)-occsvnr(ist2,jk)).gt.1d-10) then
          i=i+1
	  if (.not.req) then
            nnp(i,1,ikloc)=ist1
            nnp(i,2,ikloc)=ist2
	    if (spin_me.eq.3) then
	      nnp(i,3,ikloc)=ispn
	    else
              nnp(i,3,ikloc)=1
	    endif
            docc(i,ikloc)=occsvnr(ist1,ik)-occsvnr(ist2,jk)
	  endif
	endif
      enddo !istfv2
    enddo !istfv1
  enddo !ispn
  if (.not.req) num_nnp(ikloc)=i
  if (req) max_num_nnp=max(max_num_nnp,i)
enddo !ik

return
end

subroutine getgu(req,ngvec_me,lmaxexp,mtord,uuj,ylmgq0,sfacgq0,ngumax,ngu,gu,igu)
use modmain
implicit none
! arguments
logical, intent(in) :: req
integer, intent(in) :: ngvec_me
integer, intent(in) :: lmaxexp
integer, intent(in) :: mtord
real(8), intent(in) :: uuj(0:lmaxvr,0:lmaxvr,0:lmaxexp,mtord,mtord,natmtot,ngvec_me)
complex(8), intent(in) :: ylmgq0((lmaxexp+1)**2,ngvec_me)
complex(8), intent(in) :: sfacgq0(ngvec_me,natmtot)
integer, intent(inout) :: ngumax
integer, intent(out) :: ngu(natmtot,ngvec_me)
complex(4), intent(out) :: gu(ngumax,natmtot,ngvec_me)
integer, intent(out) :: igu(4,ngumax,natmtot,ngvec_me)

integer ig,ias,i,io1,io2,l1,m1,lm1,l2,m2,lm2,l3,m3,lm3
real(8) t1
real(8), external :: gaunt

if (req) ngumax=0 

do ig=1,ngvec_me
  do ias=1,natmtot
    i=0
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
    
return
end

subroutine getmtord(lmax,mtord)
use modmain
implicit none
integer, intent(in) :: lmax
integer, intent(out) :: mtord

integer ltmp(0:lolmax)
integer is,ilo,l,nlomaxl
 
! find maximum number of local orbitals over all l-channels
nlomaxl=0
do is=1,nspecies
  ltmp=0
  do ilo=1,nlorb(is)
    ltmp(lorbl(ilo,is))=ltmp(lorbl(ilo,is))+1
  enddo
  do l=0,lolmax
    if (l.le.lmax) then
      nlomaxl=max(nlomaxl,ltmp(l))
    endif
  enddo
enddo
mtord=apwordmax+nlomaxl

return
end


subroutine getufr(lmax,mtord,ufr)
use modmain
implicit none
! arguments
integer, intent(in) :: lmax
integer, intent(in) :: mtord
real(8), intent(out) :: ufr(nrmtmax,0:lmax,mtord,natmtot)

integer is,ia,ias,l,io,ilo
integer ordl(0:lmax)

ufr=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    ordl=0
! apw functions
    do l=0,lmax
      do io=1,apword(l,is)
        ordl(l)=ordl(l)+1
        ufr(1:nrmt(is),l,ordl(l),ias)=apwfr(1:nrmt(is),1,io,l,ias)
      enddo
    enddo
! lo functions
    do ilo=1,nlorb(is)
      l=lorbl(ilo,is)
      if (l.le.lmax) then
        ordl(l)=ordl(l)+1
        ufr(1:nrmt(is),l,ordl(l),ias)=lofr(1:nrmt(is),1,ilo,ias)
      endif
    enddo
  enddo
enddo

return
end




