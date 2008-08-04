subroutine response_chi0
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none

! number of G-vectors for matrix elements
integer ngvec_me
! q-vector in lattice coordinates
real(8) vq0l(3)
! q-vector in Cartesian coordinates
real(8) vq0c(3)
! reduced q-vector in lattice coordinates
real(8) vq0rl(3)
! reduced q-vector in Cartesian coordinates
real(8) vq0rc(3)
! index of G-vector which brings q to first BZ
integer igq0
! array for k+q indexes
integer, allocatable :: ikq(:)
! number of energy-mesh points
integer nepts
! energy mesh
complex(8), allocatable :: w(:)
! number of n,n' combinations of band indexes for each k-point
integer, allocatable :: num_nnp(:)
! maximum num_nnp over all k-points 
integer max_num_nnp
! pair of n,n' band indexes for each k-point
integer, allocatable :: nnp(:,:,:)

real(8), allocatable :: docc(:,:)
real(8), allocatable :: evalsvnr(:,:)
complex(8), allocatable :: zrhofc(:,:,:)
complex(8), allocatable :: zrhofc1(:,:)
complex(8), allocatable :: chi0_loc(:,:,:)
complex(8), allocatable :: chi0(:,:,:)
complex(8), allocatable :: mtrx1(:,:)

integer i,ik,ie,nkptnr_,ngsh_me_,i1,i2,ikloc,ig1,ig2
complex(8) wt
character*100 fname

! for parallel execution
integer, allocatable :: nkptloc(:)
integer, allocatable :: ikptloc(:,:)
integer, allocatable :: ikptiproc(:)
integer ierr,tag
integer, allocatable :: status(:)


if (iproc.eq.0) then
  write(150,'("Calculation of KS polarisability chi0")')
  write(150,*)
  write(150,'("Energy mesh parameters:")')
  write(150,'("  maximum energy [eV] : ", F7.2)')maxomega
  write(150,'("  energy step    [eV] : ", F7.2)')domega
  write(150,'("  eta            [eV] : ", F7.2)')eta
endif
  
! setup energy mesh
nepts=1+maxomega/domega
allocate(w(nepts))
do i=1,nepts
  w(i)=dcmplx(domega*(i-1),eta)/ha2ev
enddo

allocate(ikq(nkptnr))

write(fname,'("ZRHOFC[",I4.3,",",I4.3,",",I4.3,"].OUT")') &
  ivq0l(1),ivq0l(2),ivq0l(3)

if (iproc.eq.0) then
  write(150,'("Reading file ",A40)')trim(fname)
  open(160,file=trim(fname),form='unformatted',status='old')
  read(160)nkptnr_,ngsh_me_,ngvec_me,max_num_nnp,igq0
  if (nkptnr_.ne.nkptnr) then
    write(*,*)
    write(*,'("Error(response_chi0): k-mesh was changed")')
    write(*,*)
    call pstop
  endif
  if (ngsh_me_.ne.ngsh_me) then
    write(*,*)
    write(*,'("Error(response_chi0): number of G-shells for matrix elements &
      &was changed")')
    write(*,*)
    call pstop
  endif
endif
#ifdef _MPI_
call mpi_bcast(ngvec_me,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(max_num_nnp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(igq0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
allocate(num_nnp(nkptnr))
allocate(nnp(nkptnr,max_num_nnp,2))
allocate(docc(nkptnr,max_num_nnp))
if (iproc.eq.0) then
  read(160)vq0l(1:3)
  read(160)vq0rl(1:3)
  read(160)vq0c(1:3)
  read(160)vq0rc(1:3)
  do ik=1,nkptnr
    read(160)ikq(ik)
    read(160)num_nnp(ik)
    read(160)nnp(ik,1:num_nnp(ik),1:2)
    read(160)docc(ik,1:num_nnp(ik))
  enddo
endif  
#ifdef _MPI_
call mpi_bcast(vq0l,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(vq0rl,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(vq0c,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(vq0rc,3,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(num_nnp,nkptnr,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(nnp,nkptnr*max_num_nnp*2,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(docc,nkptnr*max_num_nnp,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
call mpi_bcast(ikq,nkptnr,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

! get eigen-values
allocate(evalsvnr(nstsv,nkptnr))
if (iproc.eq.0) then  
  do ik=1,nkptnr
    call getevalsv(vklnr(1,ik),evalsvnr(1,ik))
  enddo
endif
#ifdef _MPI_
call mpi_bcast(evalsvnr,nstsv*nkptnr,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif

allocate(zrhofc1(ngvec_me,max_num_nnp))
allocate(mtrx1(ngvec_me,ngvec_me))
if (iproc.eq.0) allocate(chi0(ngvec_me,ngvec_me,nepts))

! different implementations for parallel and serial execution
#ifdef _MPI_

allocate(chi0_loc(ngvec_me,ngvec_me,nepts))
allocate(status(MPI_STATUS_SIZE))
allocate(nkptloc(0:nproc-1))
allocate(ikptloc(0:nproc-1,2))
allocate(ikptiproc(nkptnr))
call splitk(nkptloc,ikptloc)
do i=0,nproc-1
  ikptiproc(ikptloc(i,1):ikptloc(i,2))=i
enddo

allocate(zrhofc(ngvec_me,max_num_nnp,nkptloc(iproc)))
do ik=1,nkptnr
! only root proc reads matrix elements from file
  if (iproc.eq.0) then
    read(160)i1,i2
    if (i1.ne.ik.or.ikq(ik).ne.i2) then
      write(*,*)
      write(*,'("Error(response_chi0): failed to read file ZRHOFC[,,].OUT")')
      write(*,*)
      call pstop
    endif
    read(160)zrhofc1(1:ngvec_me,1:num_nnp(ik))
    if (ik.le.nkptloc(0)) then
      zrhofc(:,:,ik)=zrhofc1(:,:)
    else
      tag=ik
      call mpi_send(zrhofc1,ngvec_me*max_num_nnp,MPI_DOUBLE_COMPLEX,ikptiproc(ik),tag,MPI_COMM_WORLD,ierr)
    endif
  else
    if (ik.ge.ikptloc(iproc,1).and.ik.le.ikptloc(iproc,2)) then
      ikloc=ik-ikptloc(iproc,1)+1
      tag=ik
      call mpi_recv(zrhofc1,ngvec_me*max_num_nnp,MPI_DOUBLE_COMPLEX,0,tag,MPI_COMM_WORLD,status,ierr)
      zrhofc(:,:,ikloc)=zrhofc1(:,:)
    endif
  endif
enddo
if (iproc.eq.0) then
  close(160)
  write(150,'("Finished reading matrix elements")')
endif

chi0_loc=dcmplx(0.d0,0.d0)
do ikloc=1,nkptloc(iproc)
  ik=ikptloc(iproc,1)+ikloc-1
  do i=1,num_nnp(ik)
    do ig1=1,ngvec_me
      do ig2=1,ngvec_me
        mtrx1(ig1,ig2)=zrhofc(ig1,i,ikloc)*dconjg(zrhofc(ig2,i,ikloc))
      enddo
    enddo
    do ie=1,nepts
      wt=docc(ik,i)/(evalsvnr(nnp(ik,i,1),ik)-evalsvnr(nnp(ik,i,2),ikq(ik))+w(ie))
      call zaxpy(ngvec_me**2,wt,mtrx1,1,chi0_loc(1,1,ie),1)
    enddo !ie
  enddo !i
enddo !ikloc

call mpi_reduce(chi0_loc,chi0,ngvec_me*ngvec_me*nepts,MPI_DOUBLE_COMPLEX,MPI_SUM,0,MPI_COMM_WORLD,ierr)

deallocate(chi0_loc)
deallocate(status)
deallocate(nkptloc)
deallocate(ikptloc)
deallocate(ikptiproc)
deallocate(zrhofc)

#else

chi0=dcmplx(0.d0,0.d0)
do ik=1,nkptnr
  read(160)i1,i2
  if (i1.ne.ik.or.ikq(ik).ne.i2) then
    write(*,*)
    write(*,'("Error(response_chi0): failed to read file ZRHOFC[,,].OUT")')
    write(*,*)
    call pstop
  endif
  read(160)zrhofc1(1:ngvec_me,1:num_nnp(ik))

  do i=1,num_nnp(ik)
    do ig1=1,ngvec_me
      do ig2=1,ngvec_me
        mtrx1(ig1,ig2)=zrhofc1(ig1,i)*dconjg(zrhofc1(ig2,i))
      enddo
    enddo
    do ie=1,nepts
      wt=docc(ik,i)/(evalsvnr(nnp(ik,i,1),ik)-evalsvnr(nnp(ik,i,2),ikq(ik))+w(ie))
      call zaxpy(ngvec_me**2,wt,mtrx1,1,chi0(1,1,ie),1)
    enddo !ie
  enddo !i
enddo !ik
close(160)

#endif

if (iproc.eq.0) then
  chi0=chi0/nkptnr/omega
  write(fname,'("CHI0[",I4.3,",",I4.3,",",I4.3,"].OUT")') &
    ivq0l(1),ivq0l(2),ivq0l(3)
  open(160,file=trim(fname),form='unformatted',status='replace')
  write(160)ngsh_me,ngvec_me,nepts,igq0
  write(160)w(1:nepts)
  write(160)vq0l(1:3)
  write(160)vq0rl(1:3)
  write(160)vq0c(1:3)
  write(160)vq0rc(1:3)
  do ie=1,nepts
    write(160)chi0(1:ngvec_me,1:ngvec_me,ie)
  enddo
  close(160)
endif

deallocate(w)
deallocate(ikq)
deallocate(num_nnp)
deallocate(nnp)
deallocate(docc)
deallocate(evalsvnr)
deallocate(zrhofc1)
deallocate(mtrx1)
if (iproc.eq.0) deallocate(chi0)

return
end
