subroutine wfkq(ikstep,wfsvmtloc,wfsvitloc,ngknr,igkignr,wfsvmt2, &
  wfsvit2,ngknr2,igkignr2,evecfv2,evecsv2)
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
complex(8), intent(out) :: evecfv2(nmatmax,nstfv,nspnfv)
complex(8), intent(out) :: evecsv2(nstsv,nstsv)


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
! send evecfv      
      tag=tag+1
      call mpi_isend(evecfvloc(1,1,1,ikloc),nmatmax*nstfv*nspnfv,    &
        MPI_DOUBLE_COMPLEX,rank,tag,comm_cart_100,req,ierr)
! send evecsv      
      tag=tag+1
      call mpi_isend(evecsvloc(1,1,ikloc),nstsv*nstsv,    &
        MPI_DOUBLE_COMPLEX,rank,tag,comm_cart_100,req,ierr)
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
      tag=tag+1
      call mpi_recv(evecfv2,nmatmax*nstfv*nspnfv,MPI_DOUBLE_COMPLEX, &
        rank,tag,comm_cart_100,stat,ierr)
      tag=tag+1
      call mpi_recv(evecsv2,nstsv*nstsv,MPI_DOUBLE_COMPLEX, &
        rank,tag,comm_cart_100,stat,ierr)
    endif
    if (mpi_x(1).eq.i-1.and.ip.eq.i) then
      wfsvmt2(:,:,:,:,:)=wfsvmtloc(:,:,:,:,:,ikloc)
      wfsvit2(:,:,:)=wfsvitloc(:,:,:,ikloc)
      ngknr2=ngknr(ikloc)
      igkignr2(:)=igkignr(:,ikloc)
      evecfv2(:,:,:)=evecfvloc(:,:,:,ikloc)
      evecsv2(:,:)=evecsvloc(:,:,ikloc)
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
  evecfv2(:,:,:)=evecfvloc(:,:,:,jk)
  evecsv2(:,:)=evecsvloc(:,:,jk)
#endif
return
end
