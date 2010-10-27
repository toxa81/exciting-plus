subroutine svdchi0(chi0wan,iq,w)
use modmain
use mod_addons_q

! arguments
complex(8), intent(in) :: chi0wan
integer, intent(in) :: iq
real(8), intent(in) :: w
! local variables
real(8), allocatable :: A(:,:)
real(8), allocatable :: S(:)
real(8), allocatable :: U(:,:)
real(8), allocatable :: VT(:,:)

integer i,j,iwf1,iwf2
integer lda,ldu,ldvt,lwork,info,lwmax
real(8), allocatable :: work(:)

integer nlow
real(8) t1

! do allocations
lwmax = 100000
allocate(A(nmegqwan,nmegqwan))
allocate(S(nmegqwan))
allocate(U(nmegqwan,nmegqwan))
allocate(VT(nmegqwan,nmegqwan))
allocate(work(lwmax))

! make copy of chi0wan
A = abs(chi0wan)

! lapack nonsense
lwork = lwmax
lda = nmegqwan
ldu = nmegqwan
ldvt = nmegqwan

!print *, ''
!print *, ''
!print *, ''
!print *, 'A matrix:'
!print *, '------------------------------------------------------------'
!do i=1,nmegqwan
!  print '(100F15.5)',A(i,:)
!  !write (*,"(F15.5)")(A(i,j),j=1,nmegqwan)
!enddo

call dgesvd('A','A',nmegqwan,nmegqwan,A(1,1),lda,S(1),U(1,1),ldu,VT(1,1),ldvt,work,&
  lwork,info)

nlow = -1
print *, ''
print *, ''
print *, ''
print *, 'S matrix:'
print *, '------------------------------------------------------------'
do i=1,nmegqwan
  t1 = S(1)/S(i)
  if (t1.lt.4.0) then
    nlow = i
  endif
  write(*,'(I6," value:  ",F15.5)')i,S(i)
enddo

print *, ''
print *, ''
print *, 'U matrix:'
print *, '------------------------------------------------------------'
do i=1,nlow
  print '(100F15.5)',U(:,i)
enddo
print *, ''
print *, ''
print *, ''

! deallocate
deallocate(A)
deallocate(S)
deallocate(U)
deallocate(VT)

end

