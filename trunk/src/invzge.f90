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




