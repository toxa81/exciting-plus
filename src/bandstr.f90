
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: bandstr
! !INTERFACE:
subroutine bandstr
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
integer lmax,lmmax,l,m,lm
integer ik,ispn,is,ia,ias,iv,ist,ikloc
real(8) emin,emax,sum
character(256) fname
! allocatable arrays
real(8), allocatable :: evalfv(:,:)
real(8), allocatable :: e(:,:)
! low precision for band character array saves memory
real(4), allocatable :: bc(:,:,:,:)
complex(8), allocatable :: dmat(:,:,:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
! initialise universal variables
call init0
call init1
! allocate array for storing the eigenvalues
allocate(e(nstsv,nkpt))
e=0.d0
! maximum angular momentum for band character
lmax=min(3,lmaxapw)
lmmax=(lmax+1)**2
if (task.eq.21) then
  allocate(bc(0:lmax,natmtot,nstsv,nkpt))
  bc=0.0
  allocate(dmat(lmmax,lmmax,nspinor,nspinor,nstsv))
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
end if
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
! generate muffin-tin effective magnetic fields and s.o. coupling functions
call genbeffmt
! get radial-muffint tin functions
call getufr
! get product of radial functions
call genufrp  
emin=1.d5
emax=-1.d5
allocate(evalfv(nstfv,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
! begin parallel loop over k-points
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
  write(*,'("Info(bandstr): ",I6," of ",I6," k-points")') ik,nkpt
! solve the first- and second-variational secular equations
  call seceqn(ikloc,evalfv,evecfv,evecsv)
 do ist=1,nstsv
! subtract the Fermi energy
   e(ist,ik)=evalsv(ist,ik)  !-efermi
 end do
! compute the band characters if required
  if (task.eq.21) then
! find the matching coefficients
    do ispn=1,nspnfv
      call match(ngk(ispn,ik),gkc(:,ispn,ikloc),tpgkc(:,:,ispn,ikloc), &
       sfacgk(:,:,ispn,ikloc),apwalm(:,:,:,:,ispn))
    end do
! average band character over spin and m for all atoms
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
! generate the diagonal of the density matrix
        call gendmat(.true.,.true.,0,lmax,is,ia,ngk(:,ik),apwalm,evecfv, &
         evecsv,lmmax,dmat)
        do ist=1,nstsv
          do l=0,lmax
            sum=0.d0
            do m=-l,l
              lm=idxlm(l,m)
              do ispn=1,nspinor
                sum=sum+dble(dmat(lm,lm,ispn,ispn,ist))
              end do
            end do
            bc(l,ias,ist,ik)=real(sum)
          end do
        end do
      end do
    end do
  end if
! end loop over k-points
end do
deallocate(evalfv,evecfv,evecsv)
call mpi_grid_reduce(e(1,1),nkpt*nstsv,dims=(/dim_k/))
if (task.eq.21) then
  do ik=1,nkpt
    call mpi_grid_reduce(bc(1,1,1,ik),(lmax+1)*natmtot*nstsv,dims=(/dim_k/))
  enddo
endif
do ik=1,nkpt
  do ist=1,nstsv
    emin=min(emin,e(ist,ik))
    emax=max(emax,e(ist,ik))
  enddo
enddo
emax=emax+(emax-emin)*0.5d0
emin=emin-(emax-emin)*0.5d0
if (mpi_grid_root()) then
! output the band structure
  if (task.eq.20) then
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
  else
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        write(fname,'("BAND_S",I2.2,"_A",I4.4,".OUT")') is,ia
        open(50,file=trim(fname),action='WRITE',form='FORMATTED')
        do ist=1,nstsv
          do ik=1,nkpt
! sum band character over l
            sum=0.d0
            do l=0,lmax
              sum=sum+bc(l,ias,ist,ik)
            end do
            write(50,'(2G18.10,8F12.6)') dpp1d(ik),e(ist,ik),sum, &
             (bc(l,ias,ist,ik),l=0,lmax)
          end do
          write(50,'("     ")')
        end do
        close(50)
      end do
    end do
    write(*,*)
    write(*,'("Info(bandstr):")')
    write(*,'(" band structure plot written to BAND_Sss_Aaaaa.OUT")')
    write(*,'("  for all species and atoms")')
  end if
  write(*,*)
!  write(*,'(" Fermi energy is at zero in plot")')
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
if (task.eq.21) then
  deallocate(bc)
  deallocate(apwalm)
  deallocate(dmat)
endif
return
end subroutine
!EOC
