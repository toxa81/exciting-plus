
! Copyright (C) 2002-2008 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: dos
! !INTERFACE:
subroutine dos
! !USES:
use modmain
use modtest
! !DESCRIPTION:
!   Produces a total and partial density of states (DOS) for plotting. The total
!   DOS is written to the file {\tt TDOS.OUT} while the partial DOS is written
!   to the file {\tt PDOS\_Sss\_Aaaaa.OUT} for atom {\tt aaaa} of species
!   {\tt ss}. In the case of the partial DOS, each symmetrised
!   $(l,m)$-projection is written consecutively and separated by blank lines.
!   If the global variable {\tt lmirep} is {\tt .true.}, then the density matrix
!   from which the $(l,m)$-projections are obtained is first rotated into a
!   irreducible representation basis, i.e. one that block diagonalises all the
!   site symmetry matrices in the $Y_{lm}$ basis. Eigenvalues of a quasi-random
!   matrix in the $Y_{lm}$ basis which has been symmetrised with the site
!   symmetries are written to {\tt ELMIREP.OUT}. This allows for identification
!   of the irreducible representations of the site symmetries, for example $e_g$
!   or $t_{2g}$, by the degeneracies of the eigenvalues. In the plot, spin-up is
!   made positive and spin-down negative. See the routines {\tt gendmat} and
!   {\tt brzint}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!   Parallelised and included sum over m, November 2009 (F. Cricchio)
!EOP
!BOC
implicit none
! local variables
logical tsqaz
integer lmax,lmmax
integer l0,l1,l,m,lm
integer is,ia,ias
integer nsd,ispn,jspn
integer ngknr(2),nsk(3)
integer ik,jk,igk,ist,iw,ikloc
real(8) dw,th,sps(2),vl(3),vc(3)
real(8) v1(3),v2(3),v3(3),t1
complex(8) su2(2,2),dm1(2,2),dm2(2,2)
character(256) fname
! allocatable arrays
integer, allocatable :: igkignr(:,:)
! low precision for band/spin character array saves memory
real(4), allocatable :: bc(:,:,:,:,:)
real(4), allocatable :: sc(:,:,:)
real(8), allocatable :: vgklnr(:,:,:)
real(8), allocatable :: vgkcnr(:,:,:)
real(8), allocatable :: gkcnr(:,:),tpgkcnr(:,:,:)
real(8), allocatable :: w(:),e(:,:,:),f(:,:)
real(8), allocatable :: g(:),dt(:,:),dp(:,:,:)
real(8), allocatable :: elm(:,:)
complex(8), allocatable :: ulm(:,:,:),a(:,:)
complex(8), allocatable :: dmat(:,:,:,:,:)
complex(8), allocatable :: sdmat(:,:,:)
complex(8), allocatable :: sfacgknr(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
! always use the second-variational states (avoids confusion for RDMFT)
tevecsv=.true.
! initialise universal variables
call init0
call init1
lmax=min(3,lmaxapw)
lmmax=(lmax+1)**2
if (dosssum) then
  nsd=1
else
  nsd=2
end if
if (dosmsum) then
  l0=0; l1=lmax
else
  l0=1; l1=lmmax
end if
allocate(bc(lmmax,nspinor,natmtot,nstsv,nkptnrloc))
bc=0.0
allocate(sc(nspinor,nstsv,nkptnrloc))
sc=0.0
if (lmirep) then
  allocate(elm(lmmax,natmtot))
  allocate(ulm(lmmax,lmmax,natmtot))
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
! get the eigenvalues and occupancies from file
do ik=1,nkpt
  call getevalsv(vkl(:,ik),evalsv(:,ik))
  if (dosocc) call getoccsv(vkl(:,ik),occsv(:,ik))
end do
! generate unitary matrices which convert the (l,m) basis into the irreducible
! representation basis of the symmetry group at each atomic site
if (lmirep) call genlmirep(lmax,lmmax,elm,ulm)
! compute the SU(2) operator used for rotating the density matrix to the
! desired spin-quantisation axis
v1(:)=sqados(:)
t1=sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
if (t1.le.epslat) then
  write(*,*)
  write(*,'("Error(dos): spin-quantisation axis (sqados) has zero length")')
  write(*,*)
  stop
end if
v1(:)=v1(:)/t1
if (v1(3).ge.1.d0-epslat) then
  tsqaz=.true.
else
  tsqaz=.false.
  v2(1:2)=0.d0
  v2(3)=1.d0
  call r3cross(v1,v2,v3)
! note that the spin-quantisation axis is rotated, so the density matrix should
! be rotated in the opposite direction
  th=-acos(v1(3))
  call axangsu2(v3,th,su2)
end if
allocate(igkignr(ngkmax,nspnfv))
allocate(vgklnr(3,ngkmax,nspnfv))
allocate(vgkcnr(3,ngkmax,nspnfv))
allocate(gkcnr(ngkmax,nspnfv))
allocate(tpgkcnr(2,ngkmax,nspnfv))
allocate(sfacgknr(ngkmax,natmtot,nspnfv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
allocate(dmat(lmmax,lmmax,nspinor,nspinor,nstsv))
allocate(sdmat(nspinor,nspinor,nstsv))
allocate(a(lmmax,lmmax))
! begin parallel loop over k-points
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  write(*,*)'ik=',ik
! loop over first-variational spins
  do ispn=1,nspnfv
    vl(:)=vklnr(:,ik)
    vc(:)=vkcnr(:,ik)
! spin-spiral case
    if (spinsprl) then
      if (ispn.eq.1) then
        vl(:)=vl(:)+0.5d0*vqlss(:)
        vc(:)=vc(:)+0.5d0*vqcss(:)
      else
        vl(:)=vl(:)-0.5d0*vqlss(:)
        vc(:)=vc(:)-0.5d0*vqcss(:)
      end if
    end if
! generate the G+k vectors
    call gengpvec(vl,vc,ngknr(ispn),igkignr(:,ispn),vgklnr(:,:,ispn), &
     vgkcnr(:,:,ispn),gkcnr(1,ispn),tpgkcnr(1,1,ispn))
! generate the spherical coordinates of the G+k vectors
    do igk=1,ngknr(ispn)
      call sphcrd(vgkcnr(:,igk,ispn),gkcnr(igk,ispn),tpgkcnr(:,igk,ispn))
    end do
! generate structure factors for G+k vectors
    call gensfacgp(ngknr(ispn),vgkcnr(:,:,ispn),ngkmax,sfacgknr(:,:,ispn))
! find the matching coefficients
    call match(ngknr(ispn),gkcnr(:,ispn),tpgkcnr(:,:,ispn),sfacgknr(:,:,ispn), &
     apwalm(:,:,:,:,ispn))
  end do
! get the eigenvectors from file for non-reduced k-point
  call getevecfv(vklnr(:,ik),vgklnr,evecfv)
  call getevecsv(vklnr(:,ik),evecsv)
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
! generate the density matrix
      call gendmat(.false.,.false.,0,lmax,is,ia,ngknr,apwalm,evecfv,evecsv, &
       lmmax,dmat)
! convert (l,m) part to an irreducible representation if required
      if (lmirep) then
        do ist=1,nstsv
          do ispn=1,nspinor
            do jspn=1,nspinor
              call zgemm('N','N',lmmax,lmmax,lmmax,zone,ulm(:,:,ias),lmmax, &
               dmat(:,:,ispn,jspn,ist),lmmax,zzero,a,lmmax)
              call zgemm('N','C',lmmax,lmmax,lmmax,zone,a,lmmax,ulm(:,:,ias), &
               lmmax,zzero,dmat(:,:,ispn,jspn,ist),lmmax)
            end do
          end do
        end do
      end if
! spin rotate the density matrices to desired spin-quantisation axis
      if (spinpol.and.(.not.tsqaz)) then
        do ist=1,nstsv
          do lm=1,lmmax
            dm1(:,:)=dmat(lm,lm,:,:,ist)
            call z2mm(su2,dm1,dm2)
            call z2mmct(dm2,su2,dm1)
            dmat(lm,lm,:,:,ist)=dm1(:,:)
          end do
        end do
      end if
! determine the band characters from the density matrix
      do ist=1,nstsv
        do ispn=1,nspinor
          do lm=1,lmmax
            t1=dble(dmat(lm,lm,ispn,ispn,ist))
            bc(lm,ispn,ias,ist,ikloc)=real(t1)
          end do
        end do
      end do
    end do
  end do
! compute the spin density matrices of the second-variational states
  call gensdmat(evecsv,sdmat)
! spin rotate the density matrices to desired spin-quantisation axis
  if (spinpol.and.(.not.tsqaz)) then
    do ist=1,nstsv
      call z2mm(su2,sdmat(:,:,ist),dm1)
      call z2mmct(dm1,su2,sdmat(:,:,ist))
    end do
  end if
  do ist=1,nstsv
    do ispn=1,nspinor
      t1=dble(sdmat(ispn,ispn,ist))
      sc(ispn,ist,ikloc)=real(t1)
    end do
  end do
end do
deallocate(igkignr,vgklnr,vgkcnr,gkcnr)
deallocate(tpgkcnr,sfacgknr,apwalm)
deallocate(evecfv,evecsv,dmat,sdmat,a)
allocate(w(nwdos))
allocate(e(nstsv,nkptnr,nspinor))
allocate(f(nstsv,nkptnr))
allocate(g(nwdos))
allocate(dt(nwdos,nsd))
allocate(dp(nwdos,l0:l1,nsd))
! generate energy grid
dw=(wdos(2)-wdos(1))/dble(nwdos)
do iw=1,nwdos
  w(iw)=dw*dble(iw-1)+wdos(1)
end do
! number of subdivisions used for interpolation
nsk(:)=max(ngrdos/ngridk(:),1)
! sign for spin in DOS
sps(1)=1.d0
sps(2)=-1.d0
!-------------------!
!     total DOS     !
!-------------------!
dt(:,:)=0.d0
e=0.d0
do ispn=1,nspinor
  f=0.d0
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    jk=ikmap(ivknr(1,ik),ivknr(2,ik),ivknr(3,ik))
    do ist=1,nstsv
! subtract the Fermi energy
      e(ist,ik,ispn)=evalsv(ist,jk)   !-efermi
! use diagonal of spin density matrix for weight
      f(ist,ik)=sc(ispn,ist,ikloc)
      if (dosocc) then
        f(ist,ik)=f(ist,ik)*occsv(ist,jk)
      else
        f(ist,ik)=f(ist,ik)*occmax
      end if
    end do
  end do
  call mpi_grid_reduce(f(1,1),nstsv*nkptnr,dims=(/dim_k/))
  call mpi_grid_reduce(e(1,1,ispn),nstsv*nkptnr,dims=(/dim_k/))
! integrate over the Brillouin zone
  call brzint(nsmdos,ngridk,nsk,ikmapnr,nwdos,wdos,nstsv,nstsv,e(1,1,ispn),f,g)
  if (dosssum) then
    dt(:,1)=dt(:,1)+g(:)
  else
    dt(:,ispn)=g(:)
  end if
end do
if (mpi_grid_root()) then
! output to file
  open(50,file='TDOS.OUT',action='WRITE',form='FORMATTED')
  do ispn=1,nsd
    do iw=1,nwdos
      write(50,'(2G18.10)') w(iw),dt(iw,ispn)*sps(ispn)
    end do
    write(50,'("     ")')
  end do
  close(50)
endif
!---------------------!
!     partial DOS     !
!---------------------!
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    dp(:,:,:)=0.d0
    do ispn=1,nspinor
      do l=0,lmax
        do m=-l,l
          lm=idxlm(l,m)
          f=0.d0
          do ikloc=1,nkptnrloc
            ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
            jk=ikmap(ivknr(1,ik),ivknr(2,ik),ivknr(3,ik))
            do ist=1,nstsv
              f(ist,ik)=bc(lm,ispn,ias,ist,ikloc)
              if (dosocc) then
                f(ist,ik)=f(ist,ik)*occsv(ist,jk)
              else
                f(ist,ik)=f(ist,ik)*occmax
              end if
            end do
          end do
          call mpi_grid_reduce(f(1,1),nstsv*nkptnr,dims=(/dim_k/))
          call brzint(nsmdos,ngridk,nsk,ikmapnr,nwdos,wdos,nstsv,nstsv, &
           e(:,:,ispn),f,g)
          if (dosmsum) then
            if (dosssum) then
              dp(:,l,1)=dp(:,l,1)+g(:)
            else
              dp(:,l,ispn)=dp(:,l,ispn)+g(:)
            end if
          else
            if (dosssum) then
              dp(:,lm,1)=dp(:,lm,1)+g(:)
            else
              dp(:,lm,ispn)=g(:)
            end if
          end if
! subtract from interstitial DOS
          if (dosssum) then
            dt(:,1)=dt(:,1)-g(:)
          else
            dt(:,ispn)=dt(:,ispn)-g(:)
          end if
        end do
      end do
    end do
    if (mpi_grid_root()) then
! output to file
      write(fname,'("PDOS_S",I2.2,"_A",I4.4,".OUT")') is,ia
      open(50,file=trim(fname),action='WRITE',form='FORMATTED')
      do ispn=1,nsd
        do l=l0,l1
          do iw=1,nwdos
            write(50,'(2G18.10)') w(iw),dp(iw,l,ispn)*sps(ispn)
          end do
          write(50,'("     ")')
        end do
      end do
      close(50)
    endif
  end do
end do
!------------------------------------------!
!     irreducible representations file     !
!------------------------------------------!
if (mpi_grid_root()) then
  if (lmirep) then
    open(50,file='ELMIREP.OUT',action='WRITE',form='FORMATTED')
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        write(50,*)
        write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is, &
         trim(spsymb(is)),ia
        do l=0,lmax
          do m=-l,l
            lm=idxlm(l,m)
            write(50,'(" l = ",I2,", m = ",I2,", lm= ",I3," : ",G18.10)') l,m, &
             lm,elm(lm,ias)
          end do
        end do
      end do
    end do
    close(50)
  end if
!--------------------------!
!     interstitial DOS     !
!--------------------------!
  open(50,file='IDOS.OUT',action='WRITE',form='FORMATTED')
  do ispn=1,nsd
    do iw=1,nwdos
      write(50,'(2G18.10)') w(iw),dt(iw,ispn)*sps(ispn)
    end do
    write(50,'("     ")')
  end do
  close(50)
  write(*,*)
  write(*,'("Info(dos):")')
  write(*,'(" Total density of states written to TDOS.OUT")')
  write(*,*)
  write(*,'(" Partial density of states written to PDOS_Sss_Aaaaa.OUT")')
  write(*,'(" for all species and atoms")')
  if (dosmsum) then
    write(*,'(" PDOS summed over m")')
  end if
  if (dosssum) then
    write(*,'(" PDOS summed over spin")')
  end if
  if (.not.tsqaz) then
    write(*,*)
    write(*,'(" Spin-quantisation axis : ",3G18.10)') sqados(:)
  end if
  if (lmirep) then
    write(*,*)
    write(*,'(" Eigenvalues of a random matrix in the (l,m) basis symmetrised")')
    write(*,'(" with the site symmetries written to ELMIREP.OUT for all")')
    write(*,'(" species and atoms. Degenerate eigenvalues correspond to")')
    write(*,'(" irreducible representations of each site symmetry group")')
  end if
  write(*,*)
  write(*,'(" Interstitial density of states written to IDOS.OUT")')
  write(*,*)
!  write(*,'(" Fermi energy is at zero in plot")')
!  write(*,*)
  write(*,'(" DOS units are states/Hartree/unit cell")')
  write(*,*)
endif
! write the total DOS to test file
call writetest(10,'total DOS',nv=nwdos*nsd,tol=2.d-2,rva=dt)
deallocate(bc,sc,w,e,f,g,dt,dp)
if (lmirep) deallocate(elm,ulm)
return
end subroutine
!EOC
