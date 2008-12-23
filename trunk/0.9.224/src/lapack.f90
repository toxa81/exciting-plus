      subroutine invzge(A,ndim)
      implicit none
! passed var
      integer ndim
      complex*16 A(ndim,ndim)
! local var
      integer lwork,nb,info
      real*8 ,allocatable :: work(:)
      integer ,allocatable :: ipiv(:)
      integer i,j
      integer, external :: ilaenv


      nb = ilaenv(1,'zgetri','U',ndim,-1,-1,-1)
      lwork = ndim * nb
      allocate(work(2*lwork),ipiv(ndim))

      call zgetrf(ndim,ndim,A,ndim,ipiv,info)
      if(info.ne.0) then
        write(*,*)'inverse_he_matrix: error factorization'
	stop
      endif

      call zgetri(ndim,A,ndim,ipiv,work,lwork,info)
      if(info.ne.0) then
        write(*,*)'inverse_he_matrix: error inversion'
	stop
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
  write(*,'("Warinig(invdsy) : factorization error")')
  write(*,*)
endif
call dsytri('U',n,mtrx,n,ipiv,work,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Warinig(invdsy) : inversion error")')
  write(*,*)
endif

deallocate(ipiv,work)

return
end

      subroutine diagzhe(ndim,mtrx,evalue)
      implicit   none

!---- passed var
      integer       ,intent(in   ) :: ndim
      complex*16    ,intent(inout) :: mtrx(ndim,ndim)
      real*8        ,intent(out  ) :: evalue(ndim)

!---- local var
      integer                      :: nb,lwork,inf
      real*8        ,allocatable   :: work(:),rwork(:)

      integer       ,external      :: ilaenv

      nb = ilaenv(1,'ZHETRD','U',ndim,-1,-1,-1)
      lwork = (nb+1)*ndim
      allocate(work(lwork*2))
      allocate(rwork(3*ndim+2))
      call zheev('V','U',ndim,mtrx,ndim,evalue,work,lwork,rwork,inf)
      if( inf.ne.0 ) then
        write(*,*)'DIAG_MTRX: Error finding eigenvectors.'
        stop
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
  write(*,'("Warning(diagdsy) : info = ",I4)')info
  write(*,*)
endif
deallocate(work)

return
end

subroutine isqrtzhe(ndim,mtrx,ierr)
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
  write(*,'("Error(isqrtzhe): zheev failed")')
  write(*,*)
  call pstop
endif

ierr=0
do i=1,ndim
  if (ev1(i).lt.1.d-10) then
    ierr=1
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





