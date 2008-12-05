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
