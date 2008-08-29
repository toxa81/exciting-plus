
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: dos
! !INTERFACE:
subroutine dos
! !USES:
use modmain
! !DESCRIPTION:
!   Produces a total and partial density of states (DOS) for plotting. The total
!   DOS is written to the file {\tt TDOS.OUT} while the partial DOS is written
!   to the file {\tt PDOS\_Sss\_Aaaaa.OUT} for atom {\tt aaaa} of species
!   {\tt ss}. In the case of the partial DOS, each symmetrised
!   $(l,m)$-projection is written consecutively and separated by blank lines.
!   Eigenvalues of a matrix in the $Y_{lm}$ basis which has been symmetrised
!   with the site symmetries are written to {\tt ELMSYM.OUT}. This allows for
!   identification of the irreducible representations of the site symmetries,
!   for example $e_g$ or $t_{2g}$. Spin-up is made positive and spin-down
!   negative for the partial DOS in both the collinear and non-collinear cases
!   and for the total DOS in the collinear case. See the routines {\tt bandchar}
!   and {\tt brzint}.
!
! !REVISION HISTORY:
!   Created January 2004 (JKD)
!EOP
!BOC
implicit none
! local variables
integer lmax,lmmax,l,m,lm,nsk(3),lm1,j,isym
integer ik,ispn,is,ia,ir,ias,ist,iw,i,mtord,io1,io2
real(8) t1
character(256) fname
! allocatable arrays
real(8), allocatable :: bndchr(:,:,:,:,:)
real(8), allocatable :: elmsym(:,:)
real(8), allocatable :: e(:,:,:)
real(8), allocatable :: f(:,:)
real(8), allocatable :: w(:)
real(8), allocatable :: g(:,:)
real(8), allocatable :: gp(:)
real(8), allocatable :: pdos(:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
real(8), allocatable :: ufr(:,:,:,:)
real(8), allocatable :: fr(:),gr(:),cf(:,:)
real(8), allocatable :: uu(:,:,:,:)
complex(8), allocatable :: acoeff(:,:,:,:,:)
complex(8), allocatable :: acoeff_sym(:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8) zt1(16),zt2(16)

! initialise universal variables
call init0
call init1
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
! allocate the band character array
lmax=min(3,lmaxapw)
lmmax=(lmax+1)**2
allocate(bndchr(lmmax,natmtot,nspinor,nstsv,nkpt))
allocate(elmsym(lmmax,natmtot))

allocate(fr(nrmtmax))
allocate(gr(nrmtmax))
allocate(cf(3,nrmtmax))

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

call getmtord(lmax,mtord)
allocate(ufr(nrmtmax,0:lmax,mtord,natmtot))
call getufr(lmax,mtord,ufr)
allocate(uu(0:lmax,mtord,mtord,natmtot))
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do l=0,lmax
      do io1=1,mtord
        do io2=1,mtord
          do ir=1,nrmt(is)
            fr(ir)=ufr(ir,l,io1,ias)*ufr(ir,l,io2,ias)*(spr(ir,is)**2)                                                        
          enddo
          call fderiv(-1,nrmt(is),spr(1,is),fr,gr,cf)
          uu(l,io1,io2,ias)=gr(nrmt(is))
        enddo
      enddo
    enddo 
  enddo
enddo

allocate(acoeff(lmmaxvr,mtord,natmtot,nspinor,nstsv))
allocate(acoeff_sym(lmmaxvr,mtord))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
bndchr=0.d0
do ik=1,nkpt
! get the eigenvalues/vectors from file
  call getevalsv(vkl(1,ik),evalsv(1,ik))
  call getevecfv(vkl(1,ik),vgkl(1,1,ik,1),evecfv)
  call getevecsv(vkl(1,ik),evecsv)
  call match(ngk(ik,1),gkc(1,ik,1),tpgkc(1,1,ik,1),sfacgk(1,1,ik,1),apwalm)
  call getacoeff(ngk(ik,1),mtord,apwalm,evecfv,evecsv,acoeff)
  
  do j=1,nstsv
    do ispn=1,nspinor
      do ias=1,natmtot
        acoeff_sym=dcmplx(0.d0,0.d0)
        do io1=1,mtord
          do isym=1,nsymlat
            zt2(1:16)=acoeff(1:16,io1,ias,ispn,j)
            call rotzflm(symlatc(1,1,isym),3,1,16,zt2,zt2)
            zt1=dcmplx(0.d0,0.d0)
            do lm=1,16
              do lm1=1,16
                zt1(lm)=zt1(lm)+rlm2ylm(lm1,lm)*zt2(lm1)
              enddo
            enddo
            acoeff_sym(1:16,io1)=acoeff_sym(1:16,io1)+zt1(1:16)
          enddo
        enddo
        acoeff_sym=acoeff_sym/nsymlat
        do l=0,lmax
          do m=-l,l
            lm=idxlm(l,m)
            do io1=1,mtord
              do io2=1,mtord
                bndchr(lm,ias,ispn,j,ik)=bndchr(lm,ias,ispn,j,ik) + &
                  uu(l,io1,io2,ias)*dreal(dconjg(acoeff_sym(lm,io1))*acoeff_sym(lm,io2))
              enddo
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo          
          
! compute the band character (appromximate for spin-spirals)
!  call bandchar(lmax,ik,evecfv,evecsv,lmmax,bndchr(1,1,1,1,ik),elmsym)
! compute the spin characters
  call spinchar(ik,evecsv)
end do

! generate energy grid
wdos(1)=minval(evalsv(:,:)-efermi)-0.1
wdos(2)=maxval(evalsv(:,:)-efermi)+0.1
t1=0.001d0
!t1=(wdos(2)-wdos(1))/dble(nwdos)
nwdos=1+(wdos(2)-wdos(1))/t1

! allocate local arrays
allocate(e(nstsv,nkpt,nspinor))
allocate(f(nstsv,nkpt))
allocate(w(nwdos))
allocate(g(nwdos,nspinor))
allocate(gp(nwdos))
allocate(pdos(nwdos,lmmax))

do iw=1,nwdos
  w(iw)=t1*dble(iw-1)+wdos(1)
end do
! number of subdivisions used for interpolation
do i=1,3
  nsk(i)=max(ngrdos/ngridk(i),1)
end do
! output total (spin-projected) DOS
open(50,file='TDOS.OUT',action='WRITE',form='FORMATTED')
do ispn=1,nspinor
  if (ispn.eq.1) then
    t1=1.d0
  else
    t1=-1.d0
  end if
  do ik=1,nkpt
    do ist=1,nstsv
! subtract the Fermi energy
      e(ist,ik,ispn)=evalsv(ist,ik)-efermi
! correction for scissors operator
      if (e(ist,ik,ispn).gt.0.d0) e(ist,ik,ispn)=e(ist,ik,ispn)+scissor
! use spin character for weight
      f(ist,ik)=spnchr(ispn,ist,ik)
    end do
  end do
  call brzint(nsmdos,ngridk,nsk,ikmap,nwdos,wdos,nstsv,nstsv,e(1,1,ispn),f, &
   g(1,ispn))
  do iw=1,nwdos
    write(50,'(2G18.10)') w(iw)*ha2ev,t1*g(iw,ispn)/ha2ev
  end do
  write(50,'("     ")')
end do
close(50)
! output partial DOS
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ispn=1,nspinor
      do l=0,lmax
        do m=-l,l
          lm=idxlm(l,m)
          do ik=1,nkpt
            do ist=1,nstsv
              f(ist,ik)=bndchr(lm,ias,ispn,ist,ik)
            end do
          end do
          call brzint(nsmdos,ngridk,nsk,ikmap,nwdos,wdos,nstsv,nstsv, &
           e(1,1,ispn),f,gp) 
          do iw=1,nwdos
            pdos(iw,lm) = gp(iw)/ha2ev
! interstitial DOS
            g(iw,ispn)=g(iw,ispn)-gp(iw)
          end do
        end do
      end do
      do l = 0, lmax
        write(fname,'("PDOS_S",I2.2,"_A",I4.4,"_L",I1,"_SPIN",I1,".OUT")') is,ia,l,ispn
        open(50,file=trim(fname),action='WRITE',form='FORMATTED')
        do iw = 1, nwdos
          write(50,'(8G18.10)')w(iw)*ha2ev, &
            SUM(pdos(iw,idxlm(l,-l):idxlm(l,l))),(pdos(iw,idxlm(l,m)),m=-l,l)
        enddo
        close(50)
      enddo
    end do
  end do
end do
! output eigenvalues of (l,m)-projection operator
open(50,file='ELMSYM.OUT',action='WRITE',form='FORMATTED')
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,*)
    write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is,trim(spsymb(is)),ia
    do l=0,lmax
      do m=-l,l
        lm=idxlm(l,m)
        write(50,'(" l = ",I2,", m = ",I2,", lm= ",I3," : ",G18.10)') l,m,lm, &
         elmsym(lm,ias)
      end do
    end do
  end do
end do
close(50)
! output interstitial DOS
open(50,file='IDOS.OUT',action='WRITE',form='FORMATTED')
do ispn=1,nspinor
  if (ispn.eq.1) then
    t1=1.d0
  else
    t1=-1.d0
  end if
  do iw=1,nwdos
    write(50,'(2G18.10)') w(iw)*ha2ev,t1*g(iw,ispn)/ha2ev
  end do
end do
close(50)
write(*,*)
write(*,'("Info(dos):")')
write(*,'(" total density of states written to TDOS.OUT")')
write(*,'(" partial density of states written to PDOS_Sss_Aaaaa.OUT")')
write(*,'("  for all species and atoms")')
write(*,'(" eigenvalues of a matrix in the Y_lm basis which has been")')
write(*,'("  symmetrised with the site symmetries written to ELMSYM.OUT")')
write(*,'(" interstitial density of states written to IDOS.OUT")')
write(*,*)
write(*,'(" (DOS units are states/Hartree/spin/unit cell)")')
write(*,*)
deallocate(bndchr,elmsym,e,f,w,g,gp)
deallocate(evecfv,evecsv)
return
end subroutine
!EOC
