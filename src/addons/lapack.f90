subroutine invzge(mtrx,ndim)
implicit none
! passed var
integer, intent(in) :: ndim
complex(8), intent(out) :: mtrx(ndim,ndim)
! local var
integer lwork,nb,info
real*8 ,allocatable :: work(:)
integer ,allocatable :: ipiv(:)
integer, external :: ilaenv
nb=ilaenv(1,'zgetri','U',ndim,-1,-1,-1)
lwork=ndim*nb
allocate(work(2*lwork),ipiv(ndim))
call zgetrf(ndim,ndim,mtrx,ndim,ipiv,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(invzge) zgetrf returned ",I4)')info
  call pstop
endif
call zgetri(ndim,mtrx,ndim,ipiv,work,lwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(invzge): zgetri returned ",I4)')info
  call pstop
endif
deallocate(work,ipiv)
end

subroutine invdsy(n,mtrx)
implicit none
integer, intent(in) :: n
real(8), intent(inout) :: mtrx(n,n)

real(8) t1
integer lwork,info
integer, allocatable :: ipiv(:)
real(8), allocatable :: work(:)

allocate(ipiv(n))
lwork=-1
call dsytrf('U',n,mtrx,n,ipiv,t1,lwork,info)
lwork=int(t1)+1
allocate(work(lwork))
call dsytrf('U',n,mtrx,n,ipiv,work,lwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(invdsy) dsytrf returned ",I4)')info
  call pstop
endif
call dsytri('U',n,mtrx,n,ipiv,work,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(invdsy) dsytri returned ",I4)')info
  call pstop
endif
deallocate(ipiv,work)
return
end

subroutine diagzhe(ndim,mtrx,evalue)
implicit none
integer ,intent(in) :: ndim
complex(8), intent(out) :: mtrx(ndim,ndim)
real(8), intent(out) :: evalue(ndim)
integer nb,lwork,inf
real*8, allocatable :: work(:),rwork(:)
integer, external :: ilaenv

nb=ilaenv(1,'ZHETRD','U',ndim,-1,-1,-1)
lwork=(nb+1)*ndim
allocate(work(lwork*2))
allocate(rwork(3*ndim+2))
call zheev('V','U',ndim,mtrx,ndim,evalue,work,lwork,rwork,inf)
if (inf.ne.0) then
  write(*,*)
  write(*,'("Error(diagzhe) zheev returned ",I4)')inf
  call pstop
endif
deallocate(work,rwork)
end

subroutine diagdsy(n,mtrx,eval)
implicit none
integer, intent(in) :: n
real(8), intent(inout) :: mtrx(n,n)
real(8), intent(out) :: eval(n)
integer lwork,info
real(8), allocatable :: work(:)
real(8) t1
lwork=-1
call dsyev('V','U',n,mtrx,n,eval,t1,lwork,info)
lwork=int(t1)+1
allocate(work(lwork))
call dsyev('V','U',n,mtrx,n,eval,work,lwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Warning(diagdsy) dsyev returned ",I4)')info
  write(*,*)
endif
deallocate(work)
return
end

subroutine isqrtzhe(ndim,mtrx,ierr)
use mod_mpi_grid
implicit   none
! arguments
integer, intent(in) :: ndim
integer, intent(out) :: ierr
complex(8), intent(inout) :: mtrx(ndim,ndim)

integer nb,lwork,info,i,j,n
real(8), allocatable :: work(:),rwork(:),ev(:),ev1(:)
complex(8), allocatable :: z1(:,:)

integer, external :: ilaenv

nb=ilaenv(1,'ZHETRD','U',ndim,-1,-1,-1)
lwork=(nb+1)*ndim
allocate(z1(ndim,ndim))
allocate(work(lwork*2))
allocate(rwork(3*ndim+2))
allocate(ev(ndim))
allocate(ev1(ndim))

call zheev('V','U',ndim,mtrx,ndim,ev1,work,lwork,rwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(isqrtzhe): zheev returned ",I4)')info
  write(*,*)
  call pstop
endif

ierr=0
do i=1,ndim
  if (ev1(i).lt.1.d-12) then
    ierr=i
    ev(i)=0.d0
  else
    ev(i)=1.d0/dsqrt(ev1(i))
  endif
enddo

z1(:,:)=mtrx(:,:)
mtrx=dcmplx(0.d0,0.d0)
do i=1,ndim
  do j=1,ndim
    do n=1,ndim
      mtrx(i,j)=mtrx(i,j)+z1(i,n)*dconjg(z1(j,n))*ev(n) 
    enddo
  enddo
enddo
deallocate(z1,work,rwork,ev,ev1)
return
end

subroutine isqrtdsy(n,mtrx,ierr)
use mod_mpi_grid
implicit   none
! arguments
integer, intent(in) :: n
integer, intent(out) :: ierr
real(8), intent(inout) :: mtrx(n,n)
!
integer i,j,k
real(8), allocatable :: eval(:),eval_isq(:)
real(8), allocatable :: evec(:,:)

allocate(evec(n,n))
allocate(eval(n),eval_isq(n))
call diagdsy(n,mtrx,eval)

ierr=0
do i=1,n
  if (eval(i).lt.1.d-12) then
    ierr=i
    eval_isq(i)=0.d0
  else
    eval_isq(i)=1.d0/dsqrt(eval(i))
  endif
enddo
evec(:,:)=mtrx(:,:)
mtrx=0.d0
do i=1,n
  do j=1,n
    do k=1,n
      mtrx(i,j)=mtrx(i,j)+evec(i,k)*evec(j,k)*eval_isq(k) 
    enddo
  enddo
enddo
deallocate(evec,eval,eval_isq)
return
end

subroutine diagzge(ndim,mtrx,evalue)
implicit   none
integer, intent(in) :: ndim
complex*16, intent(inout) :: mtrx(ndim,ndim)
complex*16, intent(out) :: evalue(ndim)
integer lwork,inf
complex(8) zt1
real*8, allocatable :: rwork(:)
complex(8), allocatable :: work(:)
complex(8), allocatable :: evec(:,:)
integer, external :: ilaenv

lwork=-1
call zgeev('N','V',ndim,mtrx,ndim,evalue,evec,ndim,evec,ndim,zt1,lwork,rwork,inf)
lwork=dble(zt1)+1
allocate(work(lwork))
allocate(rwork(2*ndim))
allocate(evec(ndim,ndim))
call zgeev('N','V',ndim,mtrx,ndim,evalue,evec,ndim,evec,ndim,work,lwork,rwork,inf)
if (inf.ne.0) then
  write(*,'("Error(diagzge) zgeev returned ",I4)')inf
  call pstop
endif
mtrx=evec
deallocate(work,rwork,evec)
end subroutine

!#ifdef _MAGMA_
!subroutine diagzheg(n,nv,ld,etol,a,b,eval,evec)
!implicit none
!integer, intent(in) :: n
!integer, intent(in) :: nv
!integer, intent(in) :: ld
!real(8), intent(in) :: etol
!complex(8), intent(inout) :: a(n,n)
!complex(8), intent(inout) :: b(n,n)
!real(8), intent(out) :: eval(nv)
!complex(8), intent(out) :: evec(ld,nv)
!!
!integer m,info,i,nb,lwork
!real(8) vl,vu
!integer, allocatable :: iwork(:)
!integer, allocatable :: ifail(:)
!real(8), allocatable :: w(:)
!real(8), allocatable :: rwork(:)
!complex(8), allocatable :: work(:)
!integer, external :: ilaenv 
!!
!allocate(iwork(5*n))
!allocate(ifail(n))
!allocate(w(n))
!allocate(rwork(7*n))
!allocate(work(1))
!call magmaf_zhegvx(1,'V','I','U',n,a,n,b,n,vl,vu,1,nv,etol,m,w,evec,ld, &
!  work,-1,rwork,iwork,ifail,info)
!lwork=int(dreal(work(1)))+1
!deallocate(work)
!allocate(work(lwork))
!call magmaf_zhegvx(1,'V','I','U',n,a,n,b,n,vl,vu,1,nv,etol,m,w,evec,ld, &
!  work,lwork,rwork,iwork,ifail,info)
!eval(1:nv)=w(1:nv)
!if (info.ne.0) then
!  write(*,*)
!  write(*,'("Error(diagzheg): diagonalisation failed")')
!  write(*,'("  magmaf_zhegvx returned info = ",I8)') info
!  if (info.gt.n) then
!    i=info-n
!    write(*,'(" The leading minor of the overlap matrix of order ",I8)') i
!    write(*,'("  is not positive definite")')
!    write(*,'(" Order of overlap matrix : ",I8)') n
!    write(*,*)
!  end if
!  call pstop
!end if
!deallocate(iwork,ifail,w,rwork,work)
!return
!end subroutine
!#endif


#ifdef _MAGMA_
subroutine diagzheg(n,nv,ld,etol,a,b,eval,evec)
implicit none
integer, intent(in) :: n
integer, intent(in) :: nv
integer, intent(in) :: ld
real(8), intent(in) :: etol
complex(8), intent(inout) :: a(n,n)
complex(8), intent(inout) :: b(n,n)
real(8), intent(out) :: eval(nv)
complex(8), intent(out) :: evec(ld,nv)
!
integer m,info,i,nb,lwork,lrwork,liwork
real(8) vl,vu
integer, allocatable :: iwork(:)
integer, allocatable :: isuppz(:)
real(8), allocatable :: w(:)
real(8), allocatable :: rwork(:)
complex(8), allocatable :: work(:)
integer, external :: ilaenv 
!
allocate(isuppz(2*n))
allocate(w(n))
allocate(rwork(1))
allocate(work(1))
allocate(iwork(1))
call magmaf_zhegvr(1,'V','I','U',n,a,n,b,n,vl,vu,1,nv,etol,m,w,evec,ld, &
  isuppz,work,-1,rwork,-1,iwork,-1,info)
lwork=int(dreal(work(1)))+1
lrwork=int(rwork(1))+1
liwork=iwork(1)
deallocate(rwork,iwork,work)
allocate(work(lwork))
allocate(rwork(lrwork))
allocate(iwork(liwork))
call magmaf_zhegvr(1,'V','I','U',n,a,n,b,n,vl,vu,1,nv,etol,m,w,evec,ld, &
  isuppz,work,lwork,rwork,lrwork,iwork,liwork,info)
eval(1:nv)=w(1:nv)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(diagzheg): diagonalisation failed")')
  write(*,'("  magmaf_zhegvx returned info = ",I8)') info
  if (info.gt.n) then
    i=info-n
    write(*,'(" The leading minor of the overlap matrix of order ",I8)') i
    write(*,'("  is not positive definite")')
    write(*,'(" Order of overlap matrix : ",I8)') n
    write(*,*)
  end if
  call pstop
end if
deallocate(iwork,isuppz,w,rwork,work)
return
end subroutine

#else

subroutine diagzheg(n,nv,ld,etol,a,b,eval,evec)
implicit none
integer, intent(in) :: n
integer, intent(in) :: nv
integer, intent(in) :: ld
real(8), intent(in) :: etol
complex(8), intent(inout) :: a(n,n)
complex(8), intent(inout) :: b(n,n)
real(8), intent(out) :: eval(nv)
complex(8), intent(out) :: evec(ld,nv)
!
integer m,info,i,nb,lwork
real(8) vl,vu
integer, allocatable :: iwork(:)
integer, allocatable :: ifail(:)
real(8), allocatable :: w(:)
real(8), allocatable :: rwork(:)
complex(8), allocatable :: work(:)
integer, external :: ilaenv 
!
nb=ilaenv(1,'ZHETRD','U',n,-1,-1,-1)
lwork=(nb+1)*n
allocate(iwork(5*n))
allocate(ifail(n))
allocate(w(n))
allocate(rwork(7*n))
allocate(work(lwork))
call zhegvx(1,'V','I','U',n,a,n,b,n,vl,vu,1,nv,etol,m,w,evec,ld,&
  work,lwork,rwork,iwork,ifail,info)
eval(1:nv)=w(1:nv)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(diagzheg): diagonalisation failed")')
  write(*,'("  zhegvx returned info = ",I8)') info
  if (info.gt.n) then
    i=info-n
    write(*,'(" The leading minor of the overlap matrix of order ",I8)') i
    write(*,'("  is not positive definite")')
    write(*,'(" Order of overlap matrix : ",I8)') n
    write(*,*)
  end if
  call pstop
end if
deallocate(iwork,ifail,w,rwork,work)
return
end subroutine

#endif


