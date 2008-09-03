      subroutine diag_mtrx(ndim,mtrx,evalue)
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
