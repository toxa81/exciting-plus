
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendmatlu
use modmain
implicit none
! local variables
integer ik,ispn,ist
integer is,ia,ias
real(8) t1
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: dmat(:,:,:,:,:)
integer, external :: ikglob
! allocate local arrays
allocate(evecsv(nstsv,nstsv))
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(dmat(lmmaxlu,lmmaxlu,nspinor,nspinor,nstsv))
! zero the LDA+U density matrix
dmatlu(:,:,:,:,:)=0.d0
! begin loop over k-points
do ik=1,nkptloc(iproc)
! get the eigenvectors and occupancies from file
!  call getevecfv(vkl(:,ik),vgkl(:,:,:,ik),evecfv)
!  call getevecsv(vkl(:,ik),evecsv)
!  call getoccsv(vkl(:,ik),occsv(:,ik))
! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ispn,ikglob(ik)),gkc(:,ispn,ik),tpgkc(:,:,ispn,ik), &
     sfacgk(:,:,ispn,ik),apwalm(:,:,:,:,ispn))
  end do
! begin loop over atoms and species
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      call gendmat(.false.,.false.,0,lmaxlu,is,ia,ngk(:,ikglob(ik)),apwalm, &
       evecfvloc(1,1,1,ik),evecsvloc(1,1,ik),lmmaxlu,dmat)
      do ist=1,nstsv
        t1=wkpt(ikglob(ik))*occsv(ist,ikglob(ik))
        dmatlu(:,:,:,:,ias)=dmatlu(:,:,:,:,ias)+t1*dmat(:,:,:,:,ist)
      end do
    end do
  end do
end do
call zsync(dmatlu,lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot,.true.,.true.)
! symmetrise the density matrix
call symdmat(lmaxlu,lmmaxlu,dmatlu)
deallocate(evecfv,evecsv,apwalm,dmat)
return
end subroutine

