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

real(8)                 :: vq0l(3) = (/0.d0,0.d0,0.d0/)
real(8)                 :: vq0c(3),vkq0l(3),t1
integer                 :: ik1,ik2,ist1,ist2,ik,jk,ig,ir,is
integer                 :: vgkq0l(3)
integer   , allocatable :: k1(:)
real(8)   , allocatable :: vgq0c(:,:)
real(8)   , allocatable :: gq0(:)
real(8)   , allocatable :: tpgq0(:,:)
complex(8), allocatable :: sfacgq0(:,:)
complex(8), allocatable :: ylmgq0(:,:)
real(8) ,allocatable    :: jlgq0r(:,:,:,:),jl(:)
integer ,allocatable    :: ngknr(:)
integer ,allocatable    :: igkignr(:,:)
real(8) ,allocatable    :: vgklnr(:,:,:),vgkcnr(:,:),gkcnr(:,:),tpgkcnr(:,:,:)
complex(8) ,allocatable :: sfacgknr(:,:,:)
real(8) ,external       :: r3taxi

! initialise universal variables
call init0
call init1

allocate(k1(nkptnr))
allocate(vgq0c(3,ngvec))
allocate(gq0(ngvec))
allocate(tpgq0(2,ngvec))
allocate(sfacgq0(ngvec,natmtot))
allocate(ylmgq0(lmmaxvr,ngvec))
allocate(jlgq0r(nrcmtmax,0:lmaxvr,nspecies,ngvec))
allocate(jl(0:lmaxvr))
allocate(vgklnr(3,ngkmax,nkptnr))
allocate(vgkcnr(3,ngkmax))
allocate(gkcnr(ngkmax,nkptnr))
allocate(tpgkcnr(2,ngkmax,nkptnr))
allocate(ngknr(nkptnr))
allocate(sfacgknr(ngkmax,natmtot,nkptnr))
allocate(igkignr(ngkmax,nkptnr))


!--- get q0 in Cartesian coordinates
call r3mv(bvec,vq0l,vq0c)
write(*,*)'q-vector in lattice coordinates:',vq0l
write(*,*)'q-vector in Cartesian coordinates:',vq0c

do ik = 1, nkptnr
  !--- k+q0 vector
  vkq0l(:) = vklnr(:,ik) + vq0l(:)
  vgkq0l(:) = floor(vkq0l(:))
  write(*,*)'k-point:',ik,' lattice coordinates:',vklnr(:,ik)
  write(*,*)'k+q0, lattice coordinates:',vkq0l(:)
  write(*,*)'G0, lattice coordinates:',vgkq0l(:)
  vkq0l(:) = vkq0l(:) - vgkq0l(:)
  write(*,*)'k+q0 in BZ, lattice coordinates:',vkq0l(:)
  !--- search for k' point 
  do jk = 1, nkptnr
    if (r3taxi(vklnr(:,jk),vkq0l).lt.epslat) then
      k1(ik) = jk
      goto 10
     endif
  enddo
  write(*,'("Error: k'' point not found")')
  stop
10 continue
enddo

do ik = 1, nkptnr
  write(*,*)ik,k1(ik)
enddo

!--- generate G+k vectors for entire BZ
do ik = 1, nkptnr
  call gengpvec(vklnr(:,ik),vkcnr(:,ik),ngknr(ik),igkignr,vgklnr(:,:,ik),vgkcnr,gkcnr(:,ik),tpgkcnr(:,:,ik))
  call gensfacgp(ngknr(ik),vgkcnr,ngkmax,sfacgknr(:,:,ik))
enddo


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

!--- generate Bessel functions
do ig = 1, ngvec
  do is = 1, nspecies
    do ir = 1, nrcmt(is)
      !--- |G+q0|*x
      t1 = gq0(ig)*rcmt(ir,is)
      call sbessel(lmaxvr,t1,jl)
      jlgq0r(ir,:,is,ig) = jl(:)
    enddo
  enddo
enddo



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

do ik = 1, 1 !nkptnr
  write(*,*)'k-point',ik
  
  jk = k1(ik)

  call getevecfv(vklnr(1,ik),vgklnr(:,:,ik),evecfv)
  call getevecsv(vklnr(1,ik),evecsv) 
  call match(ngknr(ik),gkcnr(:,ik),tpgkcnr(:,:,ik),sfacgknr(:,:,ik),apwalm)
  call genwfsv(.false.,ngknr(ik),igkignr(:,ik),evalsv(1,1),apwalm,evecfv, &
    evecsv,wfmt1,wfir1)
  
  call getevecfv(vklnr(1,jk),vgklnr(:,:,jk),evecfv)
  call getevecsv(vklnr(1,jk),evecsv) 
  call match(ngknr(jk),gkcnr(:,jk),tpgkcnr(:,:,jk),sfacgknr(:,:,jk),apwalm)
  call genwfsv(.false.,ngknr(jk),igkignr(:,jk),evalsv(1,1),apwalm,evecfv, &
    evecsv,wfmt2,wfir2)
  
  do ist1 = 1, 1
  do ist2 = 1, nstsv
    write(*,*)'i1=',ist1,'i2=',ist2
    call vnlrho(.true.,wfmt1(:,:,:,:,ist1),wfmt2(:,:,:,:,ist2),wfir1(:,:,ist1), &
      wfir2(:,:,ist2),zrhomt,zrhoir)
    do ig = 1, 1                                                                                                                                   
      call zrhoft(zrhomt,zrhoir,jlgq0r(:,:,:,ig),ylmgq0(:,ig),sfacgq0(ig,:))
    enddo
  enddo 
  enddo
  
enddo !ik

end
