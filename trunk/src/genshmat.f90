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

! note:
!   symmetric matrix decomposition: A=PVP^{T}, where P - matrix of eigenvectors and
!     V - diagonal matrix of eigen values v
!   from eigen-vector equation Ax=vx -> PVP^{T}x=vx -> VP^{T}x=vP^{T}x -> Vy=vy, 
!     where y=P^{T}x or x=Py; x vector in old coord.sys, y - vector in new coord.sys 
! transformation to LCS: R_m=\sum_{m'} P_{mm'} R^{loc}_{m'}, where
!   P_{mm'} matrix of eigen-vectors (stored in columns) of some matrix
!   computed in global R_m basis
!
! two more marices are required to go from R^{loc}_{lm} to Y_{lm} and back  
!
! from definition Y_{m}=\sum_{m'} rlm2ylm_{m,m'} R_{m'} in GCS:
! Y_m = \sum_{m'} rlm2ylm_{mm'}R_{m'} = \sum_{m'} rlm2ylm_{mm'} \sum_{m"} P_{m'm"} R^{loc}_{m"} = 
!    \sum_{m"} R^{loc}_{m"}  \sum_{m'} rlm2ylm_{mm'}P_{m'm"} = \sum_{m"} rlmloc2ylm_{mm"}R^{loc}_{m"}
!    where rlmloc2ylm_{mm"}=\sum_{m'} rlm2ylm_{mm'}P_{m'm"}


do ias=1,natmtot
  rlmloc2ylm(:,:,ias)=rlm2ylm(:,:)
  ylm2rlmloc(:,:,ias)=ylm2rlm(:,:)
enddo

do i=1,natlcs
  ias=iatlcs(i)
  rlmloc2ylm(:,:,ias)=dcmplx(0.d0,0.d0)
  do lm1=1,16
    do lm2=1,16
      do lm3=1,16
        rlmloc2ylm(lm1,lm2,ias)=rlmloc2ylm(lm1,lm2,ias)+rlm2ylm(lm1,lm3)*lcsrsh(lm3,lm2,i)
      enddo
    enddo
  enddo
  ylm2rlmloc(:,:,ias)=rlmloc2ylm(:,:,ias)
  call invzge(ylm2rlmloc(1,1,ias),16)
enddo

return
end

