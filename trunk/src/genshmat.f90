subroutine genshmat
use modmain
implicit none
complex(8) sqrt2,isqrt2
integer l,ias,lm1,lm2,lm3,i

ylm2rlm=dcmplx(0.d0,0.d0)
sqrt2=dcmplx(1.d0/sqrt(2.d0),0.d0)
isqrt2=dcmplx(0.d0,1.d0/sqrt(2.d0))

! construct Ylm to Rlm matrix first
!   R_{lm} = \sum_{l'm'}ylm2rlm_{lm,l'm'}Y_{l'm'}
! order of orbitals: s  y z x  xy yz 3z^2-r^2 xz x^2-y^2  f1 f2...
do l=0,3
  if (l.eq.0) then
    ylm2rlm(idxlm(0,0),idxlm(0,0))=dcmplx(1.d0,0.d0)
  else if (l.eq.1) then
    ylm2rlm(idxlm(1,-1),idxlm(1,-1))=-isqrt2
    ylm2rlm(idxlm(1,-1),idxlm(1, 1))=-isqrt2
    ylm2rlm(idxlm(1, 0),idxlm(1, 0))=dcmplx(1.d0,0.d0)
    ylm2rlm(idxlm(1, 1),idxlm(1,-1))=-sqrt2
    ylm2rlm(idxlm(1, 1),idxlm(1, 1))=sqrt2
  else if (l.eq.2) then
    ylm2rlm(idxlm(2,-2),idxlm(2,-2)) = -isqrt2
    ylm2rlm(idxlm(2,-2),idxlm(2, 2)) =  isqrt2
    ylm2rlm(idxlm(2,-1),idxlm(2,-1)) = -isqrt2
    ylm2rlm(idxlm(2,-1),idxlm(2, 1)) = -isqrt2
    ylm2rlm(idxlm(2, 0),idxlm(2, 0)) =  dcmplx(1.d0,0.d0)
    ylm2rlm(idxlm(2, 1),idxlm(2,-1)) = -sqrt2
    ylm2rlm(idxlm(2, 1),idxlm(2, 1)) =  sqrt2
    ylm2rlm(idxlm(2, 2),idxlm(2,-2)) =  sqrt2
    ylm2rlm(idxlm(2, 2),idxlm(2, 2)) =  sqrt2
  else if (l.eq.3) then
    ylm2rlm(idxlm(3,-3),idxlm(3,-3)) = -isqrt2
    ylm2rlm(idxlm(3,-3),idxlm(3, 3)) = -isqrt2
    ylm2rlm(idxlm(3,-2),idxlm(3,-2)) = -isqrt2
    ylm2rlm(idxlm(3,-2),idxlm(3, 2)) =  isqrt2
    ylm2rlm(idxlm(3,-1),idxlm(3,-1)) = -isqrt2
    ylm2rlm(idxlm(3,-1),idxlm(3, 1)) = -isqrt2
    ylm2rlm(idxlm(3, 0),idxlm(3, 0)) =  dcmplx(1.d0,0.d0)
    ylm2rlm(idxlm(3, 1),idxlm(3,-1)) = -sqrt2
    ylm2rlm(idxlm(3, 1),idxlm(3, 1)) =  sqrt2
    ylm2rlm(idxlm(3, 2),idxlm(3,-2)) =  sqrt2
    ylm2rlm(idxlm(3, 2),idxlm(3, 2)) =  sqrt2
    ylm2rlm(idxlm(3, 3),idxlm(3,-3)) = -sqrt2
    ylm2rlm(idxlm(3, 3),idxlm(3, 3)) =  sqrt2
  endif
enddo
! invert the matrix and get Rlm to Ylm transformation
rlm2ylm=ylm2rlm
call invzge(rlm2ylm,16)

! Y_{lm}=\sum_{l'm'} rlm2ylm_{lm,l'm'} R_{l'm'}
! R^{loc}_{lm}=\sum_{l'm'} D_{l'm',lm} R_{l'm'} -> R_{lm}=\sum_{l'm'} D^{-1}_{l'm',lm} R^{loc}_{l'm'}
! Y_{lm}=\sum_{l'm'} rlm2ylm_{lm,l'm'} \sum_{l"m"} D^{-1}_{l"m",l'm'} R^{loc}_{l"m"} = 
!    \sum_{l"m"} R^{loc}_{l"m"} \sum_{l'm'} rlm2ylm_{lm,l'm'} D^{-1}_{l"m",l'm'}

do ias=1,natmtot
  rlm2ylm1(:,:,ias)=rlm2ylm(:,:)
enddo

do i=1,natlcs
  ias=iatlcs(i)
  !call invdsy(16,lcsrsh(1,1,i))
  rlm2ylm1(:,:,ias)=dcmplx(0.d0,0.d0)
  do lm1=1,16
    do lm2=1,16
      do lm3=1,16
        rlm2ylm1(lm1,lm2,ias)=rlm2ylm1(lm1,lm2,ias)+rlm2ylm(lm1,lm3)*lcsrsh(lm3,lm2,i)
      enddo
    enddo
  enddo
enddo

return
end

