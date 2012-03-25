
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genapwfr
! !INTERFACE:
subroutine genapwfr1
! !USES:
use modmain
use mod_util
! !DESCRIPTION:
!   Generates the APW radial functions. This is done by integrating the scalar
!   relativistic Schr\"{o}dinger equation (or its energy deriatives) at the
!   current linearisation energies using the spherical part of the effective
!   potential. The number of radial functions at each $l$-value is given by the
!   variable {\tt apword} (at the muffin-tin boundary, the APW functions have
!   continuous derivatives up to order ${\tt apword}-1$). Within each $l$, these
!   functions are orthonormalised with the Gram-Schmidt method. The radial
!   Hamiltonian is applied to the orthonormalised functions and the results are
!   stored in the global array {\tt apwfr}.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!   Copied to equivalent atoms, February 2010 (A. Kozhevnikov and JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ja,ias,jas
integer nr,ir,nn,l,io,jo
real(8) t1
! automatic arrays
logical done(natmmax)
real(8) p1s(apwordmax)
real(8), allocatable :: vr(:),fr(:),p0(:,:),hp0(:,:)
!
allocate(vr(nrmtmax))
allocate(fr(nrmtmax))
allocate(p0(nrmtmax,apwordmax))
allocate(hp0(nrmtmax,apwordmax))

do is=1,nspecies
  nr=nrmt(is)
  done(:)=.false.
  do ia=1,natoms(is)
    if (.not.done(ia)) then
      ias=idxas(ia,is)
! use spherical part of potential
      vr(1:nr)=veffmt(1,1:nr,ias)*y00
      do l=0,lmaxapw
        do io=1,apword(l,is)
! integrate the radial Schrodinger equation
          call srrse_solve(solsc,spzn(is),apwe(io,l,ias),l,apwdm(io,l,is),nr,&
            &spr(1,is),vr,p0(1,io),hp0(1,io))
! normalise radial functions
          fr(1:nr)=p0(1:nr,io)**2
          t1=rintegrate(nr,spr(1,is),fr,m=0)
          p0(1:nr,io)=p0(1:nr,io)/sqrt(t1)
          hp0(1:nr,io)=hp0(1:nr,io)/sqrt(t1)
          p1s(io)=(p0(nr,io)-p0(nr-1,io))/(spr(nr,is)-spr(nr-1,is))
! subtract linear combination of previous vectors
          do jo=1,io-1
            fr(1:nr)=p0(1:nr,io)*p0(1:nr,jo)
            t1=rintegrate(nr,spr(1,is),fr,m=0)
            p0(1:nr,io)=p0(1:nr,io)-t1*p0(1:nr,jo)
            hp0(1:nr,io)=hp0(1:nr,io)-t1*hp0(1:nr,jo)
            p1s(io)=p1s(io)-t1*p1s(jo)
          end do
! normalise radial functions
          fr(1:nr)=p0(1:nr,io)**2
          t1=rintegrate(nr,spr(1,is),fr,m=0)
          if (t1.lt.1.d-20) then
            write(*,*)
            write(*,'("Error(genapwfr1): degenerate APW radial functions")')
            write(*,'(" for species ",I4)') is
            write(*,'(" atom ",I4)') ia
            write(*,'(" angular momentum ",I4)') l
            write(*,'(" and order ",I4)') io
            write(*,*)
            stop
          end if
          p0(1:nr,io)=p0(1:nr,io)/sqrt(t1)
          hp0(1:nr,io)=hp0(1:nr,io)/sqrt(t1)
          p1s(io)=p1s(io)/sqrt(t1)
! divide by r and store in global array
          do ir=1,nr
            t1=1.d0/spr(ir,is)
            apwfr(ir,1,io,l,ias)=t1*p0(ir,io)
            apwfr(ir,2,io,l,ias)=t1*hp0(ir,io)
          end do
! derivative at the muffin-tin surface
          apwdfr(io,l,ias)=(p1s(io)-p0(nr,io)*t1)*t1
        end do
      end do
      done(ia)=.true.
! copy to equivalent atoms
      do ja=1,natoms(is)
        if ((.not.done(ja)).and.(eqatoms(ia,ja,is))) then
          jas=idxas(ja,is)
          do l=0,lmaxapw
            do io=1,apword(l,is)
              apwfr(1:nr,1,io,l,jas)=apwfr(1:nr,1,io,l,ias)
              apwfr(1:nr,2,io,l,jas)=apwfr(1:nr,2,io,l,ias)
              apwdfr(io,l,jas)=apwdfr(io,l,ias)
            end do
          end do
          done(ja)=.true.
        end if
      end do
    end if
  end do
end do
deallocate(vr,fr,p0,hp0)
return
end subroutine
!EOC
