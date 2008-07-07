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

real(8)                 :: vq0l(3) = (/0.25d0,1.d0,0.d0/)
real(8)                 :: vq0c(3),vkq0l(3),t1
integer                 :: ik1,ik2,ist1,ist2,ik,jk,ig,ir,is,i,j,i1
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
real(8) ,allocatable    :: occsvnr(:,:)

! number of G-shells for response
integer                 :: ngsh_resp
! number of G-vectors for response (depens on ngsh_resp) 
integer                 :: ngvec_resp
! number of n,n' combinations of band-indexes for each k-point
integer ,allocatable    :: num_nnp(:)
! maximum num_nnp through all k-points 
integer                 :: max_num_nnp
! pair of n,n' band indexes for each k-point
integer ,allocatable    :: nnp(:,:,:)
! index to G-vector whcih reduces k+q0 to first BZ 
integer                 :: igkq0
! array of Fourier coefficients of complex charge density
complex(8) ,allocatable :: zrhofc(:,:)

! initialise universal variables
call init0
call init1

ngsh_resp=4

open(50,file='RESPONSE.OUT',form='formatted',status='replace')

if (task.eq.400) then
  write(50,'("Calculation of matrix elements " &
    & "<psi_{n,k}|e^{i(G+q)x}|psi_{n'',k+q}>")')
endif

! find number of G-vectors for response by given number of G-shells
ngvec_resp=1
i=1
j=1
do while (i.le.ngsh_resp)
  if (abs(gc(j+1)-gc(j)).gt.1d-10) then
    i=i+1
  endif
  j=j+1
enddo 
ngvec_resp=j-1

write(50,*)
write(50,'("Number of G-shells for response calculation:",I4)')ngsh_resp
write(50,'("Number of G-vectors for response calculation:",I4)')ngvec_resp
write(50,*)
write(50,'("G-vec.   lat.coord.    length")')
write(50,'("-----------------------------")')
do ig=1,ngvec_resp
  write(50,'(1X,I4,2X,3I4,1X,G18.10)')ig,ivg(:,ig),gc(ig)
enddo
stop


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
allocate(occsvnr(nstsv,nkptnr))

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

! check if we have enough G-shells to bring k+q0 back to first BZ
do ig = 1, ngvec_resp
  if (sum(abs(vgkq0l(:)-ivg(:,ig))).eq.0) then
    igkq0 = ig
    goto 20
  endif
enddo
write(*,*)'Not enough G-shells to bring k+q0 to first BZ'
write(*,*)'Hint: increase number of G-shells'
stop
20 continue

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

! generate structure factor for G+q0 vectors
call gensfacgp(ngvec,vgq0c,ngvec,sfacgq0)

! generate Bessel functions
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

! get occupancy of states
do ik=1,nkptnr
  call getoccsv(vklnr(1,ik),occsvnr(1,ik))
enddo

! setup n,n' stuff
! first, find the maximum size of nnp array
max_num_nnp=0
allocate(num_nnp(nkptnr))
do ik=1,nkptnr
  jk=k1(ik)
  i1=0
  do i=1,nstsv
    do j=1,nstsv
      if (abs(occsvnr(i,ik)-occsvnr(j,jk)).gt.1d-8) i1=i1+1
    enddo
  enddo
  num_nnp(ik)=i1
  max_num_nnp=max(max_num_nnp,i1)
enddo
write(*,*)"maximum number of n,n'' pairs:",max_num_nnp
allocate(nnp(nkptnr,max_num_nnp,2))
! second, setup the nnp array
do ik=1,nkptnr
  jk=k1(ik)
  i1=0
  do i=1,nstsv
    do j=1,nstsv
      if (abs(occsvnr(i,ik)-occsvnr(j,jk)).gt.1d-8) then
        i1=i1+1
        nnp(ik,i1,1)=i
        nnp(ik,i1,2)=j
      endif
    enddo
  enddo
enddo

allocate(zrhofc(ngvec_resp,max_num_nnp))


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

do ik = 1, nkptnr
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
  
  do i=1,num_nnp(ik)
    ist1=nnp(ik,i,1)
    ist2=nnp(ik,i,2)
    write(*,*)'i1=',ist1,'i2=',ist2
    call vnlrho(.true.,wfmt1(:,:,:,:,ist1),wfmt2(:,:,:,:,ist2),wfir1(:,:,ist1), &
      wfir2(:,:,ist2),zrhomt,zrhoir)
    do ig = 1, 1
      call zrhoft(zrhomt,zrhoir,jlgq0r(:,:,:,ig),ylmgq0(:,ig),sfacgq0(ig,:))
    enddo
  enddo
  
enddo !ik

end
