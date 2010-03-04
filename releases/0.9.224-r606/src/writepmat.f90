
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
integer ikloc,recl,i
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: pmat(:,:,:,:)
integer, external :: ikglob
! initialise universal variables
call init0
call init1
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
! allocate the momentum matrix elements array
allocate(pmat(3,nstsv,nstsv,nkptloc(iproc)))
! read in the density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
if (iproc.eq.0) then
  open(50,file='PMAT.OUT')
  close(50,status='DELETE')
endif
! find the record length
inquire(iolength=recl) pmat(:,:,:,1)
do ikloc=1,nkptloc(iproc)
! get the eigenvectors from file
  call getevecfv(vkl(:,ikglob(ikloc)),vgkl(:,:,:,ikloc),evecfv)
  call getevecsv(vkl(:,ikglob(ikloc)),evecsv)
! find the matching coefficients
  call match(ngk(1,ikglob(ikloc)),gkc(:,1,ikloc),tpgkc(:,:,1,ikloc), &
    sfacgk(:,:,1,ikloc),apwalm)
! calculate the momentum matrix elements
  call genpmat(ngk(1,ikglob(ikloc)),igkig(:,1,ikloc),vgkc(:,:,1,ikloc), &
    apwalm,evecfv,evecsv,pmat(1,1,1,ikloc))
end do
! write the matrix elements to direct-access file
do i=0,nproc-1
  if (i.eq.iproc) then
    open(50,file='PMAT.OUT',action='WRITE',form='UNFORMATTED',access='DIRECT', &
      recl=recl)
    do ikloc=1,nkptloc(iproc)
      write(50,rec=ikglob(ikloc)) pmat(:,:,:,ikloc)
    enddo
    close(50)
  endif
  call barrier(comm_world)
enddo
if (iproc.eq.0) then
  write(*,*)
  write(*,'("Info(writepmat):")')
  write(*,'(" momentum matrix elements written to file PMAT.OUT")')
  write(*,*)
endif
deallocate(apwalm,evecfv,evecsv,pmat)
end subroutine
!EOC

