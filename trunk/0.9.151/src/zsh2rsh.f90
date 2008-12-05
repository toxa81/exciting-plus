      subroutine zsh2rsh(zwf,rwf,ld,nr)
      use        modmain
      implicit   none
      
      integer   , intent(in   ) :: ld
      integer   , intent(in   ) :: nr
      complex(8), intent(in   ) :: zwf(ld,nrcmtmax)
      complex(8), intent(out  ) :: rwf(ld,nrcmtmax)

      complex(8)                :: sqrt2,isqrt2
      complex(8), allocatable   :: y2r(:,:)
      complex(8), allocatable   :: r2y(:,:)
      complex(8), allocatable   :: zwf_loc(:,:)
      
      integer   ,parameter      :: lmax = 3
      integer                   :: lmmax,l,m,irc,lm1,lm2
      
      lmmax = (lmax+1)**2
      allocate(y2r(lmmax,lmmax))
      allocate(r2y(lmmax,lmmax))
      allocate(zwf_loc(ld,nrcmtmax))
      
      zwf_loc = zwf
      
      if (ld.lt.lmmax) then
        write(*,*)'error: ld.lt.lmmax'
        stop
      endif

      !--- construct Ylm to Rlm matrix first
      !--- R_{lm} = \sum_{l'm'}y2r_{lm,l'm'}Y_{l'm'}
      !--- cold be done much more elegant and in more general way (using numerical technique)
      y2r = dcmplx(0.d0,0.d0)
      sqrt2  = dcmplx(1.d0/sqrt(2.d0),0.d0)
      isqrt2 = dcmplx(0.d0,1.d0/sqrt(2.d0))

      do l = 0, lmax
        if (l.eq.0) then
          y2r(idxlm(0,0),idxlm(0,0)) = dcmplx(1.d0,0.d0)
        else if (l.eq.1) then
          y2r(idxlm(1,-1),idxlm(1,-1)) = -isqrt2
          y2r(idxlm(1,-1),idxlm(1, 1)) = -isqrt2
          y2r(idxlm(1, 0),idxlm(1, 0)) =  dcmplx(1.d0,0.d0)
          y2r(idxlm(1, 1),idxlm(1,-1)) = -sqrt2
          y2r(idxlm(1, 1),idxlm(1, 1)) =  sqrt2
        else if (l.eq.2) then
          y2r(idxlm(2,-2),idxlm(2,-2)) = -isqrt2
          y2r(idxlm(2,-2),idxlm(2, 2)) =  isqrt2
          y2r(idxlm(2,-1),idxlm(2,-1)) = -isqrt2
          y2r(idxlm(2,-1),idxlm(2, 1)) = -isqrt2
          y2r(idxlm(2, 0),idxlm(2, 0)) =  dcmplx(1.d0,0.d0)
          y2r(idxlm(2, 1),idxlm(2,-1)) = -sqrt2
          y2r(idxlm(2, 1),idxlm(2, 1)) =  sqrt2
          y2r(idxlm(2, 2),idxlm(2,-2)) =  sqrt2
          y2r(idxlm(2, 2),idxlm(2, 2)) =  sqrt2
        else if (l.eq.3) then
          y2r(idxlm(3,-3),idxlm(3,-3)) = -isqrt2
          y2r(idxlm(3,-3),idxlm(3, 3)) = -isqrt2
          y2r(idxlm(3,-2),idxlm(3,-2)) = -isqrt2
          y2r(idxlm(3,-2),idxlm(3, 2)) =  isqrt2
          y2r(idxlm(3,-1),idxlm(3,-1)) = -isqrt2
          y2r(idxlm(3,-1),idxlm(3, 1)) = -isqrt2
          y2r(idxlm(3, 0),idxlm(3, 0)) =  dcmplx(1.d0,0.d0)
          y2r(idxlm(3, 1),idxlm(3,-1)) = -sqrt2
          y2r(idxlm(3, 1),idxlm(3, 1)) =  sqrt2
          y2r(idxlm(3, 2),idxlm(3,-2)) =  sqrt2
          y2r(idxlm(3, 2),idxlm(3, 2)) =  sqrt2
          y2r(idxlm(3, 3),idxlm(3,-3)) = -sqrt2
          y2r(idxlm(3, 3),idxlm(3, 3)) =  sqrt2
        else
          do m = -l, l
            y2r(idxlm(l,m),idxlm(l,m)) = dcmplx(1.d0,0.d0)
          enddo
        endif
      enddo

      !--- invert the matrix and get Rlm to Ylm transformation
      r2y = y2r
      call invzge(r2y,lmmax)
      
      rwf = dcmplx(0.d0,0.d0)
      do irc = 1, nr
        do lm1 = 1, lmmax
	  do lm2 = 1, lmmax
	    rwf(lm1,irc) = rwf(lm1,irc)+r2y(lm2,lm1)*zwf_loc(lm2,irc)
	  enddo
	enddo 
      enddo
      
      deallocate(r2y,y2r,zwf_loc)
      return 
      end
