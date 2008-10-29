subroutine genshmat
use modmain
implicit none
complex(8) sqrt2,isqrt2
integer l

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

return
end

