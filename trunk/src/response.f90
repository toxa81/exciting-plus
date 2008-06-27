subroutine response
use modmain
implicit none

! allocatable arrays                                                                                                                              
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfmt1(:,:,:,:,:),wfmt2(:,:,:,:,:)
complex(8), allocatable :: wfir1(:,:,:),wfir2(:,:,:)
complex(8), allocatable :: zrhomt(:,:,:)
complex(8), allocatable :: zrhoir(:)

integer                 :: ik1,ik2,ist1,ist2

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

ik1 = 1
ik2 = 1

call getevecfv(vkl(1,ik1),vgkl(1,1,ik1,1),evecfv)
call getevecsv(vkl(1,ik1),evecsv)
call getevalsv(vkl(1,ik1),evalsv(1,ik1))
call match(ngk(ik1,1),gkc(1,ik1,1),tpgkc(1,1,ik1,1),sfacgk(1,1,ik1,1),apwalm)                                                                         
call genwfsv(.false.,ngk(ik1,1),igkig(1,ik1,1),evalsv(1,ik1),apwalm,evecfv, &                                                                       
  evecsv,wfmt1,wfir1)

call getevecfv(vkl(1,ik2),vgkl(1,1,ik2,1),evecfv)
call getevecsv(vkl(1,ik2),evecsv)
call getevalsv(vkl(1,ik2),evalsv(1,ik2))
call match(ngk(ik2,1),gkc(1,ik2,1),tpgkc(1,1,ik2,1),sfacgk(1,1,ik2,1),apwalm)                                                                         
call genwfsv(.false.,ngk(ik2,1),igkig(1,ik2,1),evalsv(1,ik2),apwalm,evecfv, &                                                                       
  evecsv,wfmt2,wfir2)
  
do ist1 = 1, nstsv
do ist2 = 1, nstsv
  call vnlrho(.true.,wfmt1(:,:,:,:,ist1),wfmt2(:,:,:,:,ist2),wfir1(:,:,ist1), &
    wfir2(:,:,ist2),zrhomt,zrhoir)
  call zrhoft(zrhomt,zrhoir)
enddo
enddo


end
