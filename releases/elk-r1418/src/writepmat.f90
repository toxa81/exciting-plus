
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writepmat
! !INTERFACE:
subroutine writepmat
! !USES:
use modmain
! !DESCRIPTION:
!   Calculates the momentum matrix elements using routine {\tt genpmat} and
!   writes them to direct access file {\tt PMAT.OUT}.
!
! !REVISION HISTORY:
!   Created November 2003 (Sharma)
!EOP
!BOC
implicit none
! local variables
integer ik,recl,ikloc,i
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: pmat(:,:,:,:)
! initialise universal variables
call init0
call init1
! read in the density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
call getufr
if (mpi_grid_root()) then
  open(50,file='PMAT.OUT')
  close(50,status='DELETE')
endif
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(pmat(3,nstsv,nstsv,nkptloc))
! find the record length
inquire(iolength=recl) pmat(:,:,:,1)
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
! get the eigenvectors from file
  call getevecfv(vkl(:,ik),vgkl(:,:,:,ikloc),evecfv)
  call getevecsv(vkl(:,ik),evecsv)
! find the matching coefficients
  call match(ngk(1,ik),gkc(:,1,ikloc),tpgkc(:,:,1,ikloc),sfacgk(:,:,1,ikloc),&
    apwalm)
! calculate the momentum matrix elements
  call genpmat(ngk(1,ik),igkig(:,1,ikloc),vgkc(:,:,1,ikloc),apwalm,evecfv,&
    evecsv,pmat(1,1,1,ikloc))
! write the matrix elements to direct-access file
end do
! write the matrix elements to direct-access file
if (mpi_grid_side(dims=(/dim_k/))) then
  do i=0,mpi_grid_dim_size(dim_k)-1
    if (mpi_grid_dim_pos(dim_k).eq.i) then
      open(50,file='PMAT.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
        recl=recl)
      do ikloc=1,nkptloc
        ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
        write(50,rec=ik) pmat(:,:,:,ikloc)
      enddo
      close(50)
    endif
  call mpi_grid_barrier(dims=(/dim_k/))
  enddo
endif
if (mpi_grid_root()) then
  write(*,*)
  write(*,'("Info(writepmat):")')
  write(*,'(" momentum matrix elements written to file PMAT.OUT")')
  write(*,*)
endif
deallocate(evecfv,evecsv,apwalm,pmat)
end subroutine
!EOC

