
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gendmatlu
use modmain
use modldapu
implicit none
! local variables
integer ik,ispn,ist,ikloc
integer is,ia,ias
real(8) t1
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: dmat(:,:,:,:,:)
! zero the LDA+U density matrix
dmatlu(:,:,:,:,:)=0.d0
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(dmat(lmmaxlu,lmmaxlu,nspinor,nspinor,nstsv))
! begin parallel loop over k-points
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
! find the matching coefficients
  do ispn=1,nspnfv
    call match(ngk(ispn,ik),gkc(:,ispn,ikloc),tpgkc(:,:,ispn,ikloc), &
     sfacgk(:,:,ispn,ikloc),apwalm(:,:,:,:,ispn))
  end do
! begin loop over atoms and species
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      call gendmat(.false.,.false.,0,lmaxlu,is,ia,ngk(:,ik),apwalm,&
       evecfvloc(1,1,1,ikloc),evecsvloc(1,1,ikloc),lmmaxlu,dmat)
      do ist=1,nstsv
        t1=wkpt(ik)*occsv(ist,ik)
        dmatlu(:,:,:,:,ias)=dmatlu(:,:,:,:,ias)+t1*dmat(:,:,:,:,ist)
      end do
    end do
  end do
end do
deallocate(apwalm,dmat)
call mpi_grid_reduce(dmatlu(1,1,1,1,1),&
  lmmaxlu*lmmaxlu*nspinor*nspinor*natmtot,dims=(/dim_k/),all=.true.)
! symmetrise the density matrix
call symdmat(lmaxlu,lmmaxlu,dmatlu)
return
end subroutine

