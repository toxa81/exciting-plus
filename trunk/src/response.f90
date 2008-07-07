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

integer                 :: ik1,ik2,ist1,ist2,ig

real(8)                 :: vq0l(3),vq0c(3)
real(8) ,allocatable    :: vgq0c(:,:)
real(8) ,allocatable    :: gq0(:)
real(8) ,allocatable    :: tpgq0(:,:)
complex(8) ,allocatable :: sfacgq0(:,:)
integer                 :: is,ia,ias,nr,ir,l,m,lm
complex(8) ,allocatable :: ylmgq0(:,:)
real(8)                 :: t1,t2
complex(8)              :: zsum1,zsum2,zsum
real(8) ,allocatable    :: fr1(:),fr2(:),gr(:),cf(:,:)
real(8) ,allocatable    :: jlgq0r(:,:,:,:),jl(:)


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

allocate(vgq0c(3,ngvec))
allocate(gq0(ngvec))
allocate(tpgq0(2,ngvec))
allocate(sfacgq0(ngvec,natmtot))
allocate(ylmgq0(lmmaxvr,ngvec))
allocate(jlgq0r(nrcmtmax,0:lmaxvr,nspecies,ngvec))
allocate(jl(0:lmaxvr))

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

vq0l = (/0.5d0,0.5d0,0.5d0/)

!--- get q0 in lattice coordinates
call r3mv(bvec,vq0l,vq0c)

!--- generate G+q0 vectors
do ig = 1, ngvec
  vgq0c(:,ig) = vgc(:,ig) + vq0c(:)
  !--- get spherical coordinates and length of G+q0
  call sphcrd(vgq0c(:,ig),gq0(ig),tpgq0(:,ig))
  !--- generate spherical harmonics for G+q0
  call genylm(lmaxvr,tpgq0(:,ig),ylmgq0(:,ig))
enddo

!--- generate structure factor for G+q0 vectors
call gensfacgp(ngvec,vgq0c,ngvec,sfacgq0)

do ig = 1, ngvec
  do is = 1, nspecies
    nr = nrcmt(is)
    do ir = 1, nr
      !--- |G+q0|*x
      t1 = gq0(ig)*rcmt(ir,is)
      call sbessel(lmaxvr,t1,jl)
      jlgq0r(ir,:,is,ig) = jl(:)
    enddo
  enddo
enddo

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

do ist1 = 1, 50
do ist2 = 1, 100
  call vnlrho(.true.,wfmt1(:,:,:,:,ist1),wfmt2(:,:,:,:,ist2),wfir1(:,:,ist1), &
    wfir2(:,:,ist2),zrhomt,zrhoir)
  do ig = 1, 25
    call zrhoft(zrhomt,zrhoir,jlgq0r(:,:,:,ig),ylmgq0(:,ig),sfacgq0(ig,:))
  enddo
enddo
enddo


end
