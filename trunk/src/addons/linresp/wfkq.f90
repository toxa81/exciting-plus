subroutine wfkq(ikstep,wfsvmtloc,wfsvitloc,ngknr,igkignr,wfsvmt2, &
  wfsvit2,ngknr2,igkignr2,wann_c2)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
integer, intent(in) :: ikstep
complex(8), intent(in) :: wfsvmtloc(lmmaxvr,nrfmax,natmtot,nspinor,nstsv,*)
complex(8), intent(in) :: wfsvitloc(ngkmax,nspinor,nstsv,*)
integer, intent(in) :: ngknr(*)
integer, intent(in) :: igkignr(ngkmax,*)
complex(8), intent(out) :: wfsvmt2(lmmaxvr,nrfmax,natmtot,nspinor,nstsv)
complex(8), intent(out) :: wfsvit2(ngkmax,nspinor,nstsv)
integer, intent(out) :: ngknr2
integer, intent(out) :: igkignr2(ngkmax)
complex(8), intent(out) :: wann_c2(nwann,nstsv)

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
      call mpi_isend(wfsvmtloc(1,1,1,1,1,ikloc),                        &
        lmmaxvr*nrfmax*natmtot*nspinor*nstsv,                           & 
	    MPI_DOUBLE_COMPLEX,rank,tag,comm_cart_100,req,ierr)
! send interstitial part
      tag=tag+1
      call mpi_isend(wfsvitloc(1,1,1,ikloc),ngkmax*nspinor*nstsv,       & 
        MPI_DOUBLE_COMPLEX,rank,tag,comm_cart_100,req,ierr)
! send ngknr      
      tag=tag+1
      call mpi_isend(ngknr(ikloc),1,MPI_INTEGER,rank,tag,comm_cart_100, &
        req,ierr)
! send igkignr
      tag=tag+1
      call mpi_isend(igkignr(1,ikloc),ngkmax,MPI_INTEGER,rank,tag,      &
        comm_cart_100,req,ierr) 
      if (wannier) then
        tag=tag+1
        call mpi_isend(wann_c(1,1,ikloc),nwann*nstsv,                   &
          MPI_DOUBLE_COMPLEX,rank,tag,comm_cart_100,req,ierr)      
      endif
    endif
    if (mpi_x(1).eq.i-1.and.ip.ne.i) then
! source proc is ip-1
      call mpi_cart_rank(comm_cart_100,(/ip-1/),rank,ierr)
      tag=(ikstep*nproc+mpi_x(1))*10
      call mpi_recv(wfsvmt2,lmmaxvr*nrfmax*natmtot*nspinor*nstsv,       &
        MPI_DOUBLE_COMPLEX,rank,tag,comm_cart_100,stat,ierr)
      tag=tag+1
      call mpi_recv(wfsvit2,ngkmax*nspinor*nstsv,MPI_DOUBLE_COMPLEX,    &
        rank,tag,comm_cart_100,stat,ierr)
      tag=tag+1
      call mpi_recv(ngknr2,1,MPI_INTEGER,rank,tag,comm_cart_100,stat,   &
        ierr)
      tag=tag+1
      call mpi_recv(igkignr2,ngkmax,MPI_INTEGER,rank,tag,comm_cart_100, &
        stat,ierr)
      if (wannier) then
        tag=tag+1
        call mpi_recv(wann_c2,nwann*nstsv,MPI_DOUBLE_COMPLEX,rank,tag,  &
          comm_cart_100,stat,ierr)
      endif        
    endif
    if (mpi_x(1).eq.i-1.and.ip.eq.i) then
      wfsvmt2(:,:,:,:,:)=wfsvmtloc(:,:,:,:,:,ikloc)
      wfsvit2(:,:,:)=wfsvitloc(:,:,:,ikloc)
      ngknr2=ngknr(ikloc)
      igkignr2(:)=igkignr(:,ikloc)
      if (wannier) wann_c2(:,:)=wann_c(:,:,ikloc)
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
  if (wannier) wann_c2(:,:)=wann_c(:,:,jk)
#endif
return
end
