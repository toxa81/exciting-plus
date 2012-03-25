
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genlofr
! !INTERFACE:
subroutine genlofr1
! !USES:
use modmain
use mod_util
! !DESCRIPTION:
!   Generates the local-orbital radial functions. This is done by integrating
!   the scalar relativistic Schr\"{o}dinger equation (or its energy deriatives)
!   at the current linearisation energies using the spherical part of the
!   effective potential. For each local-orbital, a linear combination of
!   {\tt lorbord} radial functions is constructed such that its radial
!   derivatives up to order ${\tt lorbord}-1$ are zero at the muffin-tin radius.
!   This function is normalised and the radial Hamiltonian applied to it. The
!   results are stored in the global array {\tt lofr}.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!   Copied to equivalent atoms, February 2010 (A. Kozhevnikov and JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ja,ias,jas
integer np,nr,ir,ilo,io,jo
integer j,l,nn,info
real(8) t1
! automatic arrays
logical done(natmmax)
real(8), allocatable :: vr(:),fr(:),p0(:,:),hp0(:,:),p0s(:),hp0s(:)
! allocatable arrays
integer, allocatable :: ipiv(:)
real(8), allocatable :: xa(:),ya(:)
real(8), allocatable :: a(:,:),b(:),c(:)
! external functions
real(8) polynom
external polynom
! polynomial order
np=max(maxlorbord+1,4)
allocate(ipiv(np))
allocate(xa(np),ya(np),c(np))
allocate(a(np,np),b(np))

allocate(vr(nrmtmax))
allocate(fr(nrmtmax))
allocate(p0(nrmtmax,maxlorbord))
allocate(hp0(nrmtmax,maxlorbord))
allocate(p0s(nrmtmax))
allocate(hp0s(nrmtmax))

do is=1,nspecies
  nr=nrmt(is)
  done(:)=.false.
  do ia=1,natoms(is)
    if (.not.done(ia)) then
      ias=idxas(ia,is)
! use spherical part of potential
      vr(1:nr)=veffmt(1,1:nr,ias)*y00
      do ilo=1,nlorb(is)
        l=lorbl(ilo,is)
        do jo=1,lorbord(ilo,is)
! integrate the radial Schrodinger equation
          call srrse_solve(solsc,spzn(is),lorbe(jo,ilo,ias),l,lorbdm(jo,ilo,is),nr,&
            &spr(1,is),vr,p0(1,jo),hp0(1,jo))
          !call rschroddme(solsc,lorbdm(jo,ilo,is),l,0,lorbe(jo,ilo,ias),nprad, &
          ! nr,spr(:,is),vr,nn,p0(:,jo),p1,q0(:,jo),q1(:,jo))
! normalise radial functions
          fr(1:nr)=p0(1:nr,jo)**2
          t1=rintegrate(nr,spr(1,is),fr,m=0)
          p0(1:nr,jo)=p0(1:nr,jo)/sqrt(t1)
          hp0(1:nr,jo)=hp0(1:nr,jo)/sqrt(t1)
! set up the matrix of radial derivatives
          do j=1,np
            ir=nr-np+j
            xa(j)=spr(ir,is)
            ya(j)=p0(ir,jo)/spr(ir,is)
          end do
          do io=1,lorbord(ilo,is)
            a(io,jo)=polynom(io-1,np,xa,ya,c,rmt(is))
          end do
        end do
! set up the target vector
        b(:)=0.d0
        b(lorbord(ilo,is))=1.d0
        call dgesv(lorbord(ilo,is),1,a,np,ipiv,b,np,info)
        if (info.ne.0) then
          write(*,*)
          write(*,'("Error(genlofr1): degenerate local-orbital radial &
           &functions")')
          write(*,'(" for species ",I4)') is
          write(*,'(" atom ",I4)') ia
          write(*,'(" and local-orbital ",I4)') ilo
          write(*,'(" ZGESV returned INFO = ",I8)') info
          write(*,*)
          stop
        end if
! generate linear superposition of radial functions
        p0s(:)=0.d0
        hp0s(:)=0.d0
        do io=1,lorbord(ilo,is)
          t1=b(io)
          p0s(1:nr)=p0s(1:nr)+t1*p0(1:nr,io)
          hp0s(1:nr)=hp0s(1:nr)+t1*hp0(1:nr,io)
        end do
! normalise radial functions
        fr(1:nr)=p0s(1:nr)**2
        t1=rintegrate(nr,spr(1,is),fr,m=0)
        p0s(1:nr)=p0s(1:nr)/sqrt(t1)
        hp0s(1:nr)=hp0s(1:nr)/sqrt(t1)
! divide by r and store in global array
        do ir=1,nr
          t1=1.d0/spr(ir,is)
          lofr(ir,1,ilo,ias)=t1*p0s(ir)
          lofr(ir,2,ilo,ias)=t1*hp0s(ir)
        end do
      end do
      done(ia)=.true.
! copy to equivalent atoms
      do ja=1,natoms(is)
        if ((.not.done(ja)).and.(eqatoms(ia,ja,is))) then
          jas=idxas(ja,is)
          do ilo=1,nlorb(is)
            lofr(1:nr,1,ilo,jas)=lofr(1:nr,1,ilo,ias)
            lofr(1:nr,2,ilo,jas)=lofr(1:nr,2,ilo,ias)
          end do
          done(ja)=.true.
        end if
      end do
    end if
  end do
end do
deallocate(ipiv,xa,ya,a,b,c)
deallocate(vr,fr,p0,hp0,p0s,hp0s)
return
end subroutine
!EOC
