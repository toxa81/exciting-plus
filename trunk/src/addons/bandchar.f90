
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: bandchar
! !INTERFACE:
subroutine bandchar(dosym,lmax,ik,evecfv,evecsv,ld,bndchr)
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
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
integer, intent(in) :: ld
real(4), intent(out) :: bndchr(ld,natmtot,nspinor,nstsv)
! local variables
integer ispn,ias,io1,io2,lm
integer l,m,lm1
integer j,isym,lspl,nsym1
complex(8), allocatable :: wfsvmt(:,:,:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
! automatic arrays
complex(8) zt1(ld,nrfmax),zt2(ld),zt3
integer, external :: ikglob

allocate(wfsvmt(ld,nrfmax,natmtot,nspinor,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))

call match(ngk(1,ikglob(ik)),gkc(1,1,ik),tpgkc(1,1,1,ik),sfacgk(1,1,1,ik),apwalm)
call genwfsvmt(lmax,ld,ngk(1,ikglob(ik)),evecfv(:,:,1),evecsv,apwalm,wfsvmt)

if (dosym) then
  nsym1=nsymcrys
else
  nsym1=1
endif

bndchr=0.d0
do j=1,nstfv
  do ispn=1,nspinor
    do ias=1,natmtot
      if (dosym) then
        nsym1=nsymsite(ias)
      else
        nsym1=1
      endif
      do isym=1,nsym1
        lspl=lsplsyms(isym,ias)
        zt1=dcmplx(0.d0,0.d0)
        do io1=1,nrfmax
          call rotzflm(symlatc(1,1,lspl),3,1,ld,wfsvmt(1:ld,io1,ias,ispn,j+(ispn-1)*nstfv),zt2)
          do lm=1,ld
            do lm1=1,ld
              zt1(lm,io1)=zt1(lm,io1)+yrlm_lcs(lm1,lm,ias)*zt2(lm1)
            enddo
          enddo
        enddo !io1
        do l=0,lmax
          do m=-l,l
            lm=idxlm(l,m)
            do io1=1,nrfmax
              do io2=1,nrfmax
                bndchr(lm,ias,ispn,j+(ispn-1)*nstfv)=bndchr(lm,ias,ispn,j+(ispn-1)*nstfv) + &
                  urfprod(l,io1,io2,ias)*dreal(dconjg(zt1(lm,io1))*zt1(lm,io2))/nsym1
              enddo
            enddo
! not tested: partial contribution from APW (io2=1) or lo (io2=2)
!            zt3=zzero
!            do io1=1,nrfmax
!              zt3=zt3+dconjg(zt1(lm,io1))*urfprod(l,io1,2,ias)/nsym1
!	        enddo
!	        bndchr(lm,ias,ispn,j+(ispn-1)*nstfv)=abs(zt3)**2
          enddo !m
        enddo !l
      enddo !isym    
    enddo
  enddo
enddo

deallocate(wfsvmt,apwalm)
return
end subroutine
!EOC
