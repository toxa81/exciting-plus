      subroutine decomp_ovm(ndim,ovm,ierr)
      implicit   none

      integer       ,intent(in   ) :: ndim
      complex*16    ,intent(inout) :: ovm(ndim,ndim)
      integer       ,intent(out  ) :: ierr

      integer                      :: nb,lwork,info,i,j,n
      real*8        ,allocatable   :: work(:),rwork(:),ev(:),ev1(:)
      complex*16    ,allocatable   :: z1(:,:),z2(:,:)

      integer       ,external      :: ilaenv

      nb = ilaenv(1,'ZHETRD','U',ndim,-1,-1,-1)
      lwork = (nb+1)*ndim
      allocate(z1(ndim,ndim))
      allocate(z2(ndim,ndim))
      allocate(work(lwork*2))
      allocate(rwork(3*ndim+2))
      allocate(ev(ndim))
      allocate(ev1(ndim))

      call zheev('V','U',ndim,ovm,ndim,ev1,work,lwork,rwork,info)
      if( info.ne.0 ) then
        write(*,*)'Error(decomp_ovm): zheev failed'
	stop
      endif
!      write(lunit(1),'(255f12.6)')ev1
 
      ierr = 0
      do i = 1, ndim
        if (ev1(i).lt.1.d-10) then
	  write(*,*)'i=',i,'ev=',ev1(i)
          ierr = 2
	  ev(i) = 0.d0
	  !ev(i) = 1.d0/dsqrt(abs(ev1(i)))
	else
          ev(i) = 1.d0/dsqrt(ev1(i))
	endif
      enddo
      
!      if( ierr.eq.1 ) then
!        write(lunit(1),*)'Warning!!! ovm matrix is wrong!'
!        write(lunit(1),*)'eigen vectors of ovm:'
!	write(lunit(1),*)'real part:'
!        do i = 1, ndim
!          write(lunit(1),'(255f8.4)')(dreal(ovm(i,j)),j=1,ndim)
!        enddo
!	write(lunit(1),*)'image part:'
!        do i = 1, ndim
!          write(lunit(1),'(255f8.4)')(dimag(ovm(i,j)),j=1,ndim)
!        enddo
!        write(lunit(1),*)'eigen values of ovm:'
!        write(lunit(1),'(255f26.20)')ev1(1:ndim) 
!	write(lunit(1),*)'inverse square root of eigen values of ovm:'
!	write(lunit(1),'(255f26.20)')ev(1:ndim) 
!      endif
      
      z1(:,:) = ovm(:,:)
      do i = 1, ndim
      do j = 1, ndim
        ovm(i,j) = dcmplx(0.d0,0.d0)
	do n = 1, ndim
	  ovm(i,j) = ovm(i,j) + z1(i,n)*dconjg(z1(j,n))*ev(n)
	enddo
      enddo
      enddo

!      do i = 1, ndim
!        do j = 1, ndim
!	  z1(i,j) = ovm(i,j)*ev(j)
!	enddo
!      enddo
!
!      call zgemm('N','C',ndim,ndim,ndim,dcmplx(1.d0,0.d0),z1,          &
!                 ndim,ovm,ndim,dcmplx(0.d0,0.d0),z2,ndim)
!
!      ovm(:,:) = z2(:,:)

      deallocate(z1,z2,work,rwork,ev,ev1)

      end
