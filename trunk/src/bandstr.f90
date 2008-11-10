
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: bandstr
! !INTERFACE:
subroutine bandstr
! !USES:
use modmain
use modwann
#ifdef _MPI_
use mpi
#endif
! !DESCRIPTION:
!   Produces a band structure along the path in reciprocal-space which connects
!   the vertices in the array {\tt vvlp1d}. The band structure is obtained from
!   the second-variational eigenvalues and is written to the file {\tt BAND.OUT}
!   with the Fermi energy set to zero. If required, band structures are plotted
!   to files {\tt BAND\_Sss\_Aaaaa.OUT} for atom {\tt aaaa} of species {\tt ss},
!   which include the band characters for each $l$ component of that atom in
!   columns 4 onwards. Column 3 contains the sum over $l$ of the characters.
!   Vertex location lines are written to {\tt BANDLINES.OUT}. See routine
!   {\tt bandchar}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer lmax,lmmax,l,m,lm,ierr,i,j
integer ik,ispn,is,ia,ias,iv,ist,n
real(8) emin,emax,sum,emin0,emax0
character(256) fname
! allocatable arrays
real(8), allocatable :: evalfv(:,:)
real(4), allocatable :: bndchr(:,:,:,:,:)
real(8), allocatable :: elmsym(:,:)
real(8), allocatable :: e(:,:)
! low precision for band character array saves memory
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
integer, external :: ikglob
! initialise universal variables
call init0
call init1
if (wannier) then
  call wann_init
endif
! allocate array for storing the eigenvalues
allocate(e(nstsv,nkpt))
! maximum angular momentum for band character
lmax=min(3,lmaxapw)
lmmax=(lmax+1)**2
! allocate band character array if required
if (task.eq.21) then
  allocate(bndchr(lmmax,natmtot,nspinor,nstsv,nkpt))
end if
! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! generate the core wavefunctions and densities
call gencore
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! compute the overlap radial integrals
call olprad
! compute the Hamiltonian radial integrals
call hmlrad

if (task.eq.21.or.wannier) then
  call geturf
  call genurfprod
endif

allocate(evalfv(nstfv,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))

emin=1.d5
emax=-1.d5
e=0.d0
evalsv=0.0
bndchr=0.d0
! begin parallel loop over k-points
do ik=1,nkptloc(iproc)
  write(*,'("Info(bandstr): ",I6," of ",I6," k-points")') ikglob(ik),nkpt
! solve the first- and second-variational secular equations
  call seceqn(ik,evalfv,evecfv,evecsv)
  do ist=1,nstsv
! subtract the Fermi energy
    e(ist,ikglob(ik))=evalsv(ist,ikglob(ik))-efermi
! add scissors correction
    if (e(ist,ikglob(ik)).gt.0.d0) e(ist,ikglob(ik))=e(ist,ikglob(ik))+scissor
  end do
! compute the band characters if required
  if (task.eq.21) then
    call bandchar(.false.,lmax,ik,evecfv,evecsv,lmmax,bndchr(1,1,1,1,ikglob(ik)))
  end if
! end loop over k-points
end do
deallocate(evalfv,evecfv,evecsv)
call dsync(e,nstsv*nkpt,.true.,.false.)
if (task.eq.21) then
  do ik=1,nkpt
    call rsync(bndchr(1,1,1,1,ik),lmmax*natmtot*nspinor*nstsv,.true.,.false.)
  enddo
endif
if (wannier) then
  call dsync(wf_e,wf_dim*wann_nspins*nkpt,.true.,.false.)
endif
emin=minval(e)
emax=maxval(e)

if (iproc.eq.0) then
emax=emax+(emax-emin)*0.5d0
emin=emin-(emax-emin)*0.5d0
! output the band structure
if (task.eq.20.or.task.eq.21) then
  open(50,file='BAND.OUT',action='WRITE',form='FORMATTED')
  do ist=1,nstsv
    do ik=1,nkpt
      write(50,'(2G18.10)') dpp1d(ik),e(ist,ik)
    end do
    write(50,'("     ")')
  end do
  close(50)
  write(*,*)
  write(*,'("Info(bandstr):")')
  write(*,'(" band structure plot written to BAND.OUT")')
end if

if (wannier) then
  open(50,file='WFBAND.OUT',action='WRITE',form='FORMATTED')
  do ist=1,wf_dim
    do ik=1,nkpt
      write(50,'(2G18.10)') dpp1d(ik),wf_e(ist,1,ik)-efermi
    end do
    write(50,'("     ")')
  end do
  close(50)
endif

! output the vertex location lines
open(50,file='BANDLINES.OUT',action='WRITE',form='FORMATTED')
do iv=1,nvp1d
  write(50,'(2G18.10)') dvp1d(iv),emin
  write(50,'(2G18.10)') dvp1d(iv),emax
  write(50,'("     ")')
end do
close(50)
write(*,'(" vertex location lines written to BANDLINES.OUT")')
write(*,*)

if (task.eq.21) then
! write band-character information
  open(50,file='BANDS.OUT',action='WRITE',form='FORMATTED')
  write(50,*)lmmax,natmtot,nspinor,nstfv,nstsv,nkpt,nvp1d
    do ik = 1, nkpt
    write(50,*)dpp1d(ik)
    write(50,*)(e(ist,ik),ist=1,nstsv)
    write(50,*)((((bndchr(lm,ias,ispn,ist,ik),lm=1,lmmax), &
               ias=1,natmtot),ispn=1,nspinor),ist=1,nstsv)
  enddo
  write(50,*)wannier
  if (wannier) then
    write(50,*)wf_dim
    do ik = 1, nkpt
      write(50,*)((abs(wfc(n,i,1,ik)),n=1,wf_dim),i=1,nstfv)
    enddo
  endif
  close(50)
endif
endif

deallocate(e)
if (task.eq.21) then
  deallocate(bndchr)
endif

return
end subroutine
!EOC
