
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zrhoplot
use        modmain
implicit   none

! local variables
integer                 :: ik,ist
real(8)                 :: x,t1

integer                 :: ik1,ik2,ist1,ist2

! allocatable arrays
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfmt1(:,:,:,:,:),wfmt2(:,:,:,:,:)
complex(8), allocatable :: wfir1(:,:,:),wfir2(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)

if (nkstlist.lt.2) then
  write(*,*)
  write(*,'("Error(zrhoplot): at least two entries are required in array kstlist")')
  write(*,*)
  stop
endif
ik1  = kstlist(1,1)
ist1 = kstlist(2,1)
ik2  = kstlist(1,2)
ist2 = kstlist(2,2)

write(*,'("Plotting product of wave-function ",I3," at k-point",I3, &
  & " and wave-function ",I3," at k-pont ",I3)')ist1,ik1,ist2,ik2

! initialise universal variables
call init0
call init1

allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
allocate(wfir1(ngrtot,nspinor,nstsv))
allocate(wfir2(ngrtot,nspinor,nstsv))
allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
allocate(zrhoir(ngrtot))

write(*,*)'size of wfmt arrays: ', &
  2*lmmaxvr*nrcmtmax*natmtot*nspinor*nstsv*16.d0/1024/1024,' Mb'
write(*,*)'size of wfir arrays: ', &
  2*ngrtot*nspinor*nstsv*16.d0/1024/1024,' Mb'

! read the density and potentials from file
call readstate

! read Fermi energy from file
call readfermi

! find the new linearisation energies
call linengy

! generate the APW radial functions
call genapwfr

! generate the local-orbital radial functions
call genlofr

call getevecfv(vkl(1,ik1),vgkl(1,1,ik1,1),evecfv)
call getevecsv(vkl(1,ik1),evecsv)
call getevalsv(vkl(1,ik1),evalsv(1,ik1))
call match(ngk(ik1,1),gkc(1,ik1,1),tpgkc(1,1,ik1,1),sfacgk(1,1,ik1,1),apwalm)                                                                         
call genwfsv(.false.,ngk(ik1,1),igkig(1,ik1,1),evalsv(1,ik1),apwalm,evecfv, &                                                                       
  evecsv,wfmt1,wfir1)

if (ik1.ne.ik2) then
  call getevecfv(vkl(1,ik2),vgkl(1,1,ik2,1),evecfv)
  call getevecsv(vkl(1,ik2),evecsv)
  call getevalsv(vkl(1,ik2),evalsv(1,ik2))
  call match(ngk(ik2,1),gkc(1,ik2,1),tpgkc(1,1,ik2,1),sfacgk(1,1,ik2,1),apwalm)                                                                         
  call genwfsv(.false.,ngk(ik2,1),igkig(1,ik2,1),evalsv(1,ik2),apwalm,evecfv, &                                                                       
    evecsv,wfmt2,wfir2)
else
  wfmt2 = wfmt1
  wfir2 = wfir1
endif

!wfmt2(:,:,:,:,ist2) = dcmplx(1.d0,0.d0)
!wfir2(:,:,ist2) = dcmplx(1.d0,0.d0)

call vnlrho(.true.,wfmt1(:,:,:,:,ist1),wfmt2(:,:,:,:,ist2),wfir1(:,:,ist1), &
  wfir2(:,:,ist2),zrhomt,zrhoir)

open(50,file='ZRHO2D.OUT',action='WRITE',form='FORMATTED')
call zplot2d(50,1,lmaxvr,lmmaxvr,zrhomt,zrhoir)
close(50)

return
end subroutine

