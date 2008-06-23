
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: bandchar
! !INTERFACE:
subroutine bandchar(lmax,ik,evecfv,evecsv,ld,bndchr,elmsym)
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
integer, intent(in) :: lmax
integer, intent(in) :: ik
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
integer, intent(in) :: ld
real(8), intent(out) :: bndchr(ld,natmtot,nspinor,nstsv)
real(8), intent(out) :: elmsym(ld,natmtot)
! local variables
integer ispn,jspn,is,ia,ias,ist
integer lmmax,l,m,lm2,lm1
integer irc,i,j,n,isym,lspl
integer lwork,info
real(8) t1
complex(8) zt1,zt2
! automatic arrays
real(8) fr(nrcmtmax),gr(nrcmtmax),cf(3,nrcmtmax)
! allocatable arrays
logical, allocatable :: done(:,:)
real(8), allocatable :: rwork(:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: wfmt1(:,:,:,:)
complex(8), allocatable :: wfmt2(:,:,:)
complex(8), allocatable :: zflm(:,:)
complex(8), allocatable :: ulm(:,:,:)
complex(8), allocatable :: h0(:,:)
complex(8), allocatable :: h(:,:)
complex(8), allocatable :: work(:)
complex(8), allocatable :: y2r(:,:)
complex(8), allocatable :: r2y(:,:)
complex(8), allocatable :: wfmt3(:,:,:)
complex(8), allocatable :: wfmt4(:,:,:)
complex(8)              :: zsum
complex(8)              :: sqrt2,isqrt2

! external functions
complex(8) zdotc
external zdotc
allocate(done(nstfv,nspnfv))
lmmax=(lmax+1)**2
allocate(rwork(3*lmmax))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
if (tevecsv) allocate(wfmt1(lmmax,nrcmtmax,nstfv,nspnfv))
allocate(wfmt2(lmmax,nrcmtmax,nspinor))
allocate(zflm(lmmax,lmmax))
allocate(ulm(lmmax,lmmax,nsymlat))
allocate(h0(lmmax,lmmax))
allocate(h(lmmax,lmmax))
lwork=2*lmmax
allocate(work(lwork))
allocate(y2r(lmmax,lmmax))
allocate(r2y(lmmax,lmmax))
allocate(wfmt3(lmmax,nrcmtmax,nspinor))
allocate(wfmt4(lmmax,nrcmtmax,nspinor))

bcsym = .false.

!--- construct Ylm to Rlm matrix first
!--- R_{lm} = \sum_{l'm'}y2r_{lm,l'm'}Y_{l'm'}
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
      lm1=idxlm(l,m)
      y2r(lm1,lm1) = dcmplx(1.d0,0.d0)
    enddo
  endif
enddo

!--- invert the matrix and get Rlm to Ylm transformation
r2y = y2r
call invzge(r2y,lmmax)

! find the matching coefficients
do ispn=1,nspnfv
  call match(ngk(ik,ispn),gkc(1,ik,ispn),tpgkc(1,1,ik,ispn), &
   sfacgk(1,1,ik,ispn),apwalm(1,1,1,1,ispn))
end do

if (bcsym) then
! construct (l,m) rotation matrix for each lattice symmetry
  zflm(:,:)=0.d0
  do lm1=1,lmmax
    zflm(lm1,lm1)=1.d0
  end do
  do isym=1,nsymlat
    call rotzflm(symlatc(1,1,isym),lmax,lmmax,lmmax,zflm,ulm(1,1,isym))
  end do
! set up quasi-random symmetric matrix h0
  n=0
  do lm1=1,lmmax
    do lm2=lm1,lmmax
      n=n+1
      h0(lm1,lm2)=dble(n)
      h0(lm2,lm1)=h0(lm1,lm2)
    end do
  end do
else
! set h to the unit matrix if no symmetrisation is required
  h(:,:)=0.d0
  do lm1=1,lmmax
    h(lm1,lm1)=1.d0
  end do
  elmsym(1:lmmax,:)=0.d0
end if

bndchr = 0.d0

! begin loop over species
do is=1,nspecies
! begin loop over atoms
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    if (bcsym) then
! symmetrise h0 with site symmetries
      h(:,:)=0.d0
      do isym=1,nsymsite(ias)
! spatial rotation element in lattice point group
        lspl=lsplsyms(isym,ias)
! apply (l,m) symmetry matrix as H -> H + U*H0*conjg(U')
        call zgemm('N','N',lmmax,lmmax,lmmax,zone,ulm(1,1,lspl),lmmax,h0, &
         lmmax,zzero,zflm,lmmax)
        call zgemm('N','C',lmmax,lmmax,lmmax,zone,zflm,lmmax,ulm(1,1,lspl), &
         lmmax,zone,h,lmmax)
      end do
! block diagonalise h
      do l=0,lmax
        n=2*l+1
        lm1=idxlm(l,-l)
        call zheev('V','U',n,h(lm1,lm1),lmmax,elmsym(lm1,ias),work,lwork, &
         rwork,info)
      end do
    end if
    
    done(:,:)=.false.
    do j=1,nstsv
      if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
        wfmt2(:,:,:)=0.d0
        i=0
        do ispn=1,nspinor
          if (spinsprl) then
            jspn=ispn
          else
            jspn=1
          end if
          do ist=1,nstfv
            i=i+1
            zt1=evecsv(i,j)
            if (abs(dble(zt1))+abs(aimag(zt1)).gt.epsocc) then
              if (.not.done(ist,jspn)) then
                call wavefmt(lradstp,lmax,is,ia,ngk(ik,jspn), &
                 apwalm(1,1,1,1,jspn),evecfv(1,ist,jspn),lmmax, &
                 wfmt1(1,1,ist,jspn))
                done(ist,jspn)=.true.
              end if
! add to spinor wavefunction
              n=lmmax*nrcmt(is)
              call zaxpy(n,zt1,wfmt1(1,1,ist,jspn),1,wfmt2(1,1,ispn),1)
            end if
          end do
        end do
      else
! spin-unpolarised wavefunction
        call wavefmt(lradstp,lmax,is,ia,ngk(ik,1),apwalm,evecfv(1,j,1),lmmax, &
         wfmt2)
      end if
      
      do isym=1,nsymlat
        
        do ispn = 1, nspinor
          do irc = 1, nrcmt(is)
            call rotzflm(symlatc(1,1,isym),lmax,1,lmmax,wfmt2(:,irc,ispn),wfmt4(:,irc,ispn))
          enddo
        enddo
      
        wfmt3 = dcmplx(0.d0,0.d0)
        do ispn = 1, nspinor
          do irc = 1, nrcmt(is)
	    do lm1 = 1, lmmax
	      do lm2 = 1, lmmax
	        wfmt3(lm1,irc,ispn) = wfmt3(lm1,irc,ispn)+r2y(lm2,lm1)*wfmt4(lm2,irc,ispn)
	      enddo
	    enddo 
	  enddo
        enddo
      
! determine the band characters
        do ispn=1,nspinor
          do l=0,lmax
            n=2*l+1
            lm1=idxlm(l,-l)
            do m=-l,l
              lm2=idxlm(l,m)
              do irc=1,nrcmt(is)
! project wavefunction onto eigenvectors of h
                zt1=zdotc(n,h(lm1,lm2),1,wfmt3(lm1,irc,ispn),1)
                t1=dble(zt1)**2+aimag(zt1)**2
                fr(irc)=t1*rcmt(irc,is)**2
              end do
! perform radial integral
              call fderiv(-1,nrcmt(is),rcmt(1,is),fr,gr,cf)
              bndchr(lm2,ias,ispn,j)=bndchr(lm2,ias,ispn,j) + gr(nrcmt(is))
            end do
          end do
        end do
      
      end do !nsymlat
    
    end do !j
! end loop over atoms
  end do
! end loop over species
end do

bndchr = bndchr/nsymlat 

deallocate(done,rwork,apwalm,wfmt2)
if (tevecsv) deallocate(wfmt1)
deallocate(zflm,ulm,h0,h,work,wfmt3,wfmt4,r2y)
return
end subroutine
!EOC

