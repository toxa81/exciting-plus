subroutine getmaxgnt(lmaxexp,maxgnt)
use modmain
implicit none
integer, intent(in) :: lmaxexp
integer, intent(out) :: maxgnt
integer l1,l2,l3,m1,m2,m3
integer nrf1,nrf2
real(8), external :: gaunt
real(8) t1
! estimate the maximum number of Gaunt-like coefficients 
maxgnt=0
do l1=0,lmaxapw
  nrf1=maxval(nufr(l1,:))
  do l2=0,lmaxapw
    nrf2=maxval(nufr(l2,:))
    do m1=-l1,l1
      do m2=-l2,l2
        t1=0.d0
        do l3=0,lmaxexp
          do m3=-l3,l3
            t1=t1+abs(gaunt(l2,l1,l3,m2,m1,m3))
          enddo
        enddo
        if (t1.gt.1d-16) maxgnt=maxgnt+nrf1*nrf2
      enddo
    enddo
  enddo
enddo
return
end
      