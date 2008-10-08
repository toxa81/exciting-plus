
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: bandchar
! !INTERFACE:
subroutine bandchar(dosym,lmax,ik,mtord,evecfv,evecsv,ld,bndchr,uu)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   lmax   : maximum angular momentum (in,integer)
!   ik     : k-point number (in,integer)
!   evecfv : first-variational eigenvectors (in,complex(nmatmax,nstfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
!   ld     : leading dimension (in,integer)
!   bndchr : band character (out,real(ld,natmtot,nspinor,nstsv))
!   elmsym : eigenvalues of a matrix in the Y_lm basis which has been
!            symmetrised with the site symmetries (out,real(ld,natmtot))
! !DESCRIPTION:
!   Returns the so-called ``band characters'' of the second-variational states.
!   These are given by
!   $$ \xi^{i{\bf p}\alpha}_{lm\sigma}=\int^{R_{\alpha}}_0\left|\sum_j
!    C^i_{j\sigma}\sum_{m'}v^{\alpha*}_{lmm'}\Psi^{j{\bf p}\alpha}_{lm'}(r)
!    \right|^2r^2dr $$
!   where $\Psi^{j{\bf p}\alpha}_{lm'}$ are the $r$-dependent $(l,m)$-components
!   of the first-variational muffin-tin wavefunction for state $j$,
!   $k$-point ${\bf p}$ and atom $\alpha$; and $C_{j\sigma}^{i}$ are the
!   coefficients of spinor component $\sigma$ of the $i$th second-variational
!   state. In order to obtain a physically relevant $m$ projection, the vector
!   $v^{\alpha}_{lmm'}$ is taken to be the $m$th eigenvector of a random
!   Hermitian matrix of dimension $2l+1$ which has been symmetrised with site
!   symmetries $S^{\alpha}$ in the spherical harmonic basis:
!   $$ h=\sum_i S^{\alpha}_i h_0(S^{\alpha}_i)^{-1}. $$
!   Thus the degeneracy of the eigenvalues of $h$ will determine the irreducible
!   representations of the site symmetry group. These eigenvalues are returned
!   in the array {\tt elmsym}. If the global variable {\tt bcsym} is
!   {\tt .false.} then the band characters refer to non-symmetrised
!   $m$-projections, i.e. $v^{\alpha}_{lmm'}=\delta_{mm'}$. Band characters give
!   an indication of the spin, atomistic and $(l,m)$-strengths of each state.
!   See routines {\tt seceqnsv} and {\tt wavefmt}.
!
! !REVISION HISTORY:
!   Created December 2003 (JKD)
!   Fixed problem with second-variational states, November 2006 (JKD)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: dosym
integer, intent(in) :: lmax
integer, intent(in) :: ik
integer, intent(in) :: mtord
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
integer, intent(in) :: ld
real(4), intent(out) :: bndchr(ld,natmtot,nspinor,nstsv)
real(8), intent(in) :: uu(0:lmax,mtord,mtord,natmtot)
! local variables
integer ispn,jspn,is,ia,ias,ist,io1,io2,lm,ikglob
integer l,m,lm2,lm1
integer irc,i,j,n,isym,lspl,nsym1
integer lwork,info
real(8) t1
complex(8), allocatable :: acoeff(:,:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
! automatic arrays
real(8) fr(nrcmtmax),gr(nrcmtmax),cf(3,nrcmtmax)
complex(8) zt1(ld,mtord),zt2(ld)

ikglob=ikptloc(iproc,1)+ik-1

allocate(acoeff(ld,mtord,natmtot,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))

call match(ngk(ikglob,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1),apwalm)
call getacoeff(lmax,ld,ngk(ikglob,1),mtord,apwalm,evecfv,evecsv,acoeff)

if (dosym) then
  nsym1=nsymcrys
else
  nsym1=1
endif

bndchr=0.d0
do j=1,nstfv
  do ispn=1,nspinor
    do ias=1,natmtot
      do isym=1,nsym1
        zt1=dcmplx(0.d0,0.d0)
        do io1=1,mtord
          call rotzflm(symlatc(1,1,lsplsymc(isym)),3,1,ld,acoeff(1:ld,io1,ias,j+(ispn-1)*nstfv),zt2)
          do lm=1,ld
            do lm1=1,ld
              zt1(lm,io1)=zt1(lm,io1)+rlm2ylm(lm1,lm)*zt2(lm1)
            enddo
          enddo
        enddo !io1
        do l=0,lmax
          do m=-l,l
            lm=idxlm(l,m)
            do io1=1,mtord
              do io2=1,mtord
                bndchr(lm,ias,ispn,j)=bndchr(lm,ias,ispn,j) + &
                  uu(l,io1,io2,ias)*dreal(dconjg(zt1(lm,io1))*zt1(lm,io2))
              enddo
            enddo
          enddo
        enddo
      enddo !isym    
    enddo
  enddo
enddo
bndchr=bndchr/nsym1

deallocate(acoeff,apwalm)
return
end subroutine
!EOC







subroutine getry
use modmain
implicit none
complex(8) sqrt2,isqrt2
integer l

ylm2rlm=dcmplx(0.d0,0.d0)
sqrt2=dcmplx(1.d0/sqrt(2.d0),0.d0)
isqrt2=dcmplx(0.d0,1.d0/sqrt(2.d0))

! construct Ylm to Rlm matrix first
!   R_{lm} = \sum_{l'm'}ylm2rlm_{lm,l'm'}Y_{l'm'}
! order of orbitals: s  y z x  xy xz 3z^2-r^2 yz x^2-y^2  f1 f2...
do l=0,3
  if (l.eq.0) then
    ylm2rlm(idxlm(0,0),idxlm(0,0)) = dcmplx(1.d0,0.d0)
  else if (l.eq.1) then
    ylm2rlm(idxlm(1,-1),idxlm(1,-1)) = -isqrt2
    ylm2rlm(idxlm(1,-1),idxlm(1, 1)) = -isqrt2
    ylm2rlm(idxlm(1, 0),idxlm(1, 0)) =  dcmplx(1.d0,0.d0)
    ylm2rlm(idxlm(1, 1),idxlm(1,-1)) = -sqrt2
    ylm2rlm(idxlm(1, 1),idxlm(1, 1)) =  sqrt2
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
