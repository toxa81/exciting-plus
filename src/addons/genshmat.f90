subroutine genshmat
use modmain
implicit none
complex(8) sqrt2,isqrt2
integer l,ias,lm1,lm2,lm3,i,m1,m2

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
! transformation to LPS: R_m=\sum_{m'} P_{mm'} R^{loc}_{m'}, where
!   P_{mm'} matrix of eigen-vectors (stored in columns) of some matrix
!   computed in global R_m basis
!
! two more marices are required to go from R^{loc}_{lm} to Y_{lm} and back  
!
! from definition Y_{m}=\sum_{m'} yrlm_{m,m'} R_{m'} in GCS:
! Y_m = \sum_{m'} yrlm_{mm'}R_{m'} = \sum_{m'} yrlm_{mm'} \sum_{m"} P_{m'm"} R^{loc}_{m"} = 
!    \sum_{m"} R^{loc}_{m"}  \sum_{m'} yrlm_{mm'}P_{m'm"} = \sum_{m"} yrlm_lps{mm"}R^{loc}_{m"}
!    where yrlm_lps{mm"}=\sum_{m'} yrlm_{mm'}P_{m'm"}


do ias=1,natmtot
  yrlm_lps(:,:,ias)=yrlm(:,:)
  rylm_lps(:,:,ias)=rylm(:,:)
enddo

do i=1,natlps
  ias=iatlps(i)
  yrlm_lps(:,:,ias)=dcmplx(0.d0,0.d0)
  do lm1=1,16
    do lm2=1,16
      do lm3=1,16
        yrlm_lps(lm1,lm2,ias)=yrlm_lps(lm1,lm2,ias)+yrlm(lm1,lm3)*lpsrsh(lm3,lm2,i)
      enddo
    enddo
  enddo
  rylm_lps(:,:,ias)=yrlm_lps(:,:,ias)
  call invzge(rylm_lps(1,1,ias),16)
enddo

! new implementation
if (allocated(zdsht)) deallocate(zdsht)
allocate(zdsht(lmmaxapw,lmmaxapw))
zdsht=zzero
if (allocated(dzsht)) deallocate(dzsht)
allocate(dzsht(lmmaxapw,lmmaxapw))
dzsht=zzero
! Mathematica code to generate real spherical harmonics
!b[m1_, m2_] := 
! If[m1 == 0, 1, 
!  If[m1 < 0 && m2 < 0, -I/Sqrt[2], 
!   If[m1 < 0 && m2 > 0, (-1)^m2*I/Sqrt[2], 
!    If[m1 > 0 && m2 < 0, (-1)^m1/Sqrt[2], 
!     If[m1 > 0 && m2 > 0, 1/Sqrt[2]]]]]]
!a[m1_, m2_] := If[Abs[m1] == Abs[m2], b[m1, m2], 0]
!R[l_, m_, t_, p_] := 
! Sum[a[m, m1]*SphericalHarmonicY[l, m1, t, p], {m1, -l, l}]
do l=0,lmaxapw
  do m1=-l,l
    do m2=-l,l
      if (abs(m1).eq.abs(m2)) then
        if (m1.eq.0) zdsht(idxlm(l,m1),idxlm(l,m1))=zone
        if (m1.lt.0.and.m2.lt.0) zdsht(idxlm(l,m1),idxlm(l,m2))=-zi/sqrt(2.d0)
        if (m1.lt.0.and.m2.gt.0) zdsht(idxlm(l,m1),idxlm(l,m2))=(-1)**m2*zi/sqrt(2.d0)        
        if (m1.gt.0.and.m2.lt.0) zdsht(idxlm(l,m1),idxlm(l,m2))=(-1)**m1/sqrt(2.d0)        
        if (m1.gt.0.and.m2.gt.0) zdsht(idxlm(l,m1),idxlm(l,m2))=1.d0/sqrt(2.d0)        
      endif
    enddo
  enddo
enddo
! R_{L}=\sum_{L'} M_{L,L'} Y_{L'} so Y_{L}=\sum_{L'} (M^{-1})_{L,L'} R_{L'}
dzsht=zdsht
call invzge(dzsht,lmmaxapw)
return
end

