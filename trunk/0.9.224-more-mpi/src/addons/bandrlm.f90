
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: bandstr
! !INTERFACE:
subroutine bandrlm
! !USES:
use modmain
! !DESCRIPTION:
!   Produces a band structure along the path in reciprocal-space which connects
!   the vertices in the array {\tt vvlp1d}. The band structure is obtained from
!   the second-variational eigenvalues and is written to the file {\tt BAND.OUT}
!   with the Fermi energy set to zero. If required, band structures are plotted
!   to files {\tt BAND\_Sss\_Aaaaa.OUT} for atom {\tt aaaa} of species {\tt ss},
!   which include the band characters for each $l$ component of that atom in
!   columns 4 onwards. Column 3 contains the sum over $l$ of the characters.
!   Vertex location lines are written to {\tt BANDLINES.OUT}.
!
! !REVISION HISTORY:
!   Created June 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
integer lmax,lmmax,l,m,lm,i,j,n
integer ik,ispn,is,ia,ias,iv,ist
real(8) emin,emax,sum
character(256) fname
! allocatable arrays
real(8), allocatable :: evalfv(:,:)
real(8), allocatable :: e(:,:)
! low precision for band character array saves memory
real(4), allocatable :: bc(:,:,:,:,:)
complex(8), allocatable :: dmat(:,:,:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
integer, external :: ikglob
! initialise universal variables
call init0
call init1
! allocate array for storing the eigenvalues
allocate(e(nstsv,nkpt))
! maximum angular momentum for band character
lmax=min(3,lmaxapw)
lmmax=(lmax+1)**2
allocate(bc(lmmax,natmtot,nspinor,nstsv,nkpt))
allocate(evalfv(nstfv,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
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
call geturf
call genurfprod
! begin parallel loop over k-points
e=0.d0
bc=0.d0
do ik=1,nkptloc(iproc)
  write(*,'("Info(bandstr): ",I6," of ",I6," k-points")') ikglob(ik),nkpt
! solve the first- and second-variational secular equations
  call seceqn(ik,evalfv,evecfv,evecsv)
  do ist=1,nstsv
! subtract the Fermi energy
    e(ist,ikglob(ik))=evalsv(ist,ikglob(ik)) !-efermi
  end do
! compute the band characters if required
  call bandchar(.false.,lmax,ik,evecfv,evecsv,lmmax,bc(1,1,1,1,ikglob(ik)))
! end loop over k-points
end do
deallocate(evalfv,evecfv,evecsv) 
call dsync(e,nstsv*nkpt,.true.,.false.)
if (wannier) call dsync(wann_e,wann_nmax*wann_nspin*nkpt,.true.,.false.)
do ik=1,nkpt
  call rsync(bc(1,1,1,1,ik),lmmax*natmtot*nspinor*nstsv,.true.,.false.)
  call barrier
enddo
emin=minval(e)
emax=maxval(e)
emax=emax+(emax-emin)*0.5d0
emin=emin-(emax-emin)*0.5d0
if (iproc.eq.0) then
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
  if (wannier) then
    open(50,file='WANN_BAND.OUT',action='WRITE',form='FORMATTED')
    do ispn=1,wann_nspin
      write(50,'("# spin : ",I1)')ispn
      do ist=1,nwann(ispn)
        do ik=1,nkpt
          write(50,'(2G18.10)') dpp1d(ik),wann_e(ist,ispn,ik)
        end do
        write(50,'("     ")')
      end do
    enddo !ispn
  endif
  open(50,file='BNDCHR.OUT',action='WRITE',form='FORMATTED')
  write(50,*)lmmax,natmtot,nspinor,nstfv,nstsv,nkpt,nvp1d
  write(50,*)efermi
  do ik = 1, nkpt
    write(50,*)dpp1d(ik)
    write(50,*)(e(ist,ik),ist=1,nstsv)
    write(50,*)((((bc(lm,ias,ispn,ist,ik),lm=1,lmmax), &
               ias=1,natmtot),ispn=1,nspinor),ist=1,nstsv)
  enddo
  write(50,*)wannier
  if (wannier) then
    write(50,*)wann_nmax
  endif
endif

if (wannier) then
  do i=0,nproc-1
    if (i.eq.iproc) then
      open(50,file='BNDCHR.OUT',action='WRITE',form='FORMATTED',status='OLD',position='APPEND')
      do ik=1,nkptloc(iproc)
        write(50,*)((abs(wann_c(n,j,1,ik)),n=1,nwann(1)),j=1,nstfv)
      enddo
      close(50)
    endif
    call barrier
  enddo
endif

if (iproc.eq.0) then
write(*,*)
write(*,'(" Fermi energy is at zero in plot")')
! output the vertex location lines
open(50,file='BANDLINES.OUT',action='WRITE',form='FORMATTED')
do iv=1,nvp1d
  write(50,'(2G18.10)') dvp1d(iv),emin
  write(50,'(2G18.10)') dvp1d(iv),emax
  write(50,'("     ")')
end do
close(50)
write(*,*)
write(*,'(" vertex location lines written to BANDLINES.OUT")')
write(*,*)
endif
deallocate(e)
if (task.eq.21) deallocate(bc)
return
end subroutine
!EOC
