subroutine genshmat
use modmain
implicit none
complex(8) sqrt2,isqrt2
integer l,ias,lm1,lm2,lm3,i

rylm=dcmplx(0.d0,0.d0)
sqrt2=dcmplx(1.d0/sqrt(2.d0),0.d0)
isqrt2=dcmplx(0.d0,1.d0/sqrt(2.d0))

! construct Ylm to Rlm matrix first
!   R_{lm} = \sum_{l'm'}rylm_{lm,l'm'}Y_{l'm'}
! order of orbitals: s  y z x  xy yz 3z^2-r^2 xz x^2-y^2  f1 f2...
do l=0,3
  if (l.eq.0) then
    rylm(idxlm(0,0),idxlm(0,0))=dcmplx(1.d0,0.d0)
  else if (l.eq.1) then
    rylm(idxlm(1,-1),idxlm(1,-1))=-isqrt2
    rylm(idxlm(1,-1),idxlm(1, 1))=-isqrt2
    rylm(idxlm(1, 0),idxlm(1, 0))=dcmplx(1.d0,0.d0)
    rylm(idxlm(1, 1),idxlm(1,-1))=-sqrt2
    rylm(idxlm(1, 1),idxlm(1, 1))=sqrt2
  else if (l.eq.2) then
    rylm(idxlm(2,-2),idxlm(2,-2)) = -isqrt2
    rylm(idxlm(2,-2),idxlm(2, 2)) =  isqrt2
    rylm(idxlm(2,-1),idxlm(2,-1)) = -isqrt2
    rylm(idxlm(2,-1),idxlm(2, 1)) = -isqrt2
    rylm(idxlm(2, 0),idxlm(2, 0)) =  dcmplx(1.d0,0.d0)
    rylm(idxlm(2, 1),idxlm(2,-1)) = -sqrt2
    rylm(idxlm(2, 1),idxlm(2, 1)) =  sqrt2
    rylm(idxlm(2, 2),idxlm(2,-2)) =  sqrt2
    rylm(idxlm(2, 2),idxlm(2, 2)) =  sqrt2
  else if (l.eq.3) then
    rylm(idxlm(3,-3),idxlm(3,-3)) = -isqrt2
    rylm(idxlm(3,-3),idxlm(3, 3)) = -isqrt2
    rylm(idxlm(3,-2),idxlm(3,-2)) = -isqrt2
    rylm(idxlm(3,-2),idxlm(3, 2)) =  isqrt2
    rylm(idxlm(3,-1),idxlm(3,-1)) = -isqrt2
    rylm(idxlm(3,-1),idxlm(3, 1)) = -isqrt2
    rylm(idxlm(3, 0),idxlm(3, 0)) =  dcmplx(1.d0,0.d0)
    rylm(idxlm(3, 1),idxlm(3,-1)) = -sqrt2
    rylm(idxlm(3, 1),idxlm(3, 1)) =  sqrt2
    rylm(idxlm(3, 2),idxlm(3,-2)) =  sqrt2
    rylm(idxlm(3, 2),idxlm(3, 2)) =  sqrt2
    rylm(idxlm(3, 3),idxlm(3,-3)) = -sqrt2
    rylm(idxlm(3, 3),idxlm(3, 3)) =  sqrt2
  endif
enddo
! invert the matrix and get Rlm to Ylm transformation
yrlm=rylm
call invzge(yrlm,16)

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
! from definition Y_{m}=\sum_{m'} yrlm_{m,m'} R_{m'} in GCS:
! Y_m = \sum_{m'} yrlm_{mm'}R_{m'} = \sum_{m'} yrlm_{mm'} \sum_{m"} P_{m'm"} R^{loc}_{m"} = 
!    \sum_{m"} R^{loc}_{m"}  \sum_{m'} yrlm_{mm'}P_{m'm"} = \sum_{m"} yrlm_lcs{mm"}R^{loc}_{m"}
!    where yrlm_lcs{mm"}=\sum_{m'} yrlm_{mm'}P_{m'm"}


do ias=1,natmtot
  yrlm_lcs(:,:,ias)=yrlm(:,:)
  rylm_lcs(:,:,ias)=rylm(:,:)
enddo

do i=1,natlcs
  ias=iatlcs(i)
  yrlm_lcs(:,:,ias)=dcmplx(0.d0,0.d0)
  do lm1=1,16
    do lm2=1,16
      do lm3=1,16
        yrlm_lcs(lm1,lm2,ias)=yrlm_lcs(lm1,lm2,ias)+yrlm(lm1,lm3)*lcsrsh(lm3,lm2,i)
      enddo
    enddo
  enddo
  rylm_lcs(:,:,ias)=yrlm_lcs(:,:,ias)
  call invzge(rylm_lcs(1,1,ias),16)
enddo

return
end

