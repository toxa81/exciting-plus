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

real(8)                 :: vq0l(3) = (/0.5d0,0.5d0,0.5d0/)
real(8)                 :: vq0c(3),vkq0l(3),t1,vq0rl(3),vq0rc(3)
integer                 :: ik1,ik2,ist1,ist2,ik,jk,ig,ir,is,i,j,i1,i2,i3,ig1,ie
integer                 :: vgkq0l(3),vgq0l(3)
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
real                    :: cpu0,cpu1
complex(8) ,allocatable :: chi0(:,:,:)
real(8) ,allocatable    :: evalsvnr(:,:)
complex(8) ,allocatable :: w(:)
! number of G-shells for response
integer                 :: ngsh_resp
! number of G-vectors for response (depens on ngsh_resp) 
integer                 :: ngvec_resp
! number of n,n' combinations of band-indexes for each k-point
integer ,allocatable    :: num_nnp(:)
! maximum num_nnp over all k-points 
integer                 :: max_num_nnp
! pair of n,n' band indexes for each k-point
integer ,allocatable    :: nnp(:,:,:)
! index to G-vector whcih reduces q0 to first BZ 
integer                 :: igq0
! array of Fourier coefficients of complex charge density
complex(8) ,allocatable :: zrhofc(:,:)

integer                 :: nepts
real(8)                 :: emax

real(8) ,allocatable    :: docc(:,:)
integer    ,allocatable :: gshell(:)

! initialise universal variables
call init0
call init1

ngsh_resp=1

open(50,file='RESPONSE.OUT',form='formatted',status='replace')

if (task.eq.400) then
  write(50,'("Calculation of matrix elements <psi_{n,k}|e^{i(G+q)x}|psi_{n'',k+q}>")')

! find number of G-vectors by given number of G-shells
allocate(gshell(ngvec))
ngvec_resp=1
i=1
j=1
do while (i.le.ngsh_resp)
  gshell(j)=i
  if (abs(gc(j+1)-gc(j)).gt.epslat) then
    i=i+1
  endif
  j=j+1
enddo 
ngvec_resp=j-1

write(50,*)
write(50,'("Number of G-shells for response calculation  :",I4)')ngsh_resp
write(50,'("Number of G-vectors for response calculation :",I4)')ngvec_resp
write(50,*)
write(50,'("  G-vec.       lat.coord.      length(1/a.u.) shell")')
write(50,'(" ---------------------------------------------------")')
do ig=1,ngvec_resp
  write(50,'(2X,I4,4X,3I5,4X,F12.6,5x,I4)')ig,ivg(:,ig),gc(ig),gshell(ig)
enddo

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

! find G-vector which brings q0 to first BZ
vgq0l(:)=floor(vq0l(:))

! reduce q0 vector fo first BZ
vq0rl(:)=vq0l(:)-vgq0l(:)

! check if we have enough G-shells to bring q0 back to first BZ
do ig=1,ngvec_resp
  if (sum(abs(vgq0l(:)-ivg(:,ig))).eq.0) then
    igq0=ig
    goto 20
  endif
enddo
write(*,*)
write(*,'("Error(response): not enough G-vectors to reduce q-vector to first BZ")')
write(*,'(" Increase number of G-shells")')
write(*,*)
stop
20 continue

! get q0 and reduced q0 in Cartesian coordinates
call r3mv(bvec,vq0l,vq0c)
call r3mv(bvec,vq0rl,vq0rc)

write(50,*)
write(50,'("q-vector in lattice coordinates              : ",3G18.10)')vq0l
write(50,'("q-vector in Cartesian coordinates (1/a.u.)   : ",3G18.10)')vq0c
write(50,'("q-vector length (1/a.u.)                     : ",G18.10)')sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)
write(50,'("q-vector length (1/A)                        : ",G18.10)')sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)/0.529177d0
write(50,'("G-vector to reduce q to first BZ (lat.coord.): ",3I4)')vgq0l
write(50,'("Index of G-vector                            : ",I4)')igq0
write(50,'("Reduced q-vector (lat.coord.)                : ",3G18.10)')vq0rl

! find k+q and reduce them to first BZ (this is required to utilize the 
!   periodical property of Bloch-states: |k>=|k+K>, where K is
!   any vector of reciprocal lattice)
do ik=1,nkptnr
! k+q vector
  vkq0l(:)=vklnr(:,ik)+vq0l(:)
! K vector
  vgkq0l(:)=floor(vkq0l(:))
!-  write(*,*)'k-point:',ik,' lattice coordinates:',vklnr(:,ik)
!-  write(*,*)'k+q0, lattice coordinates:',vkq0l(:)
!-  write(*,*)'G0, lattice coordinates:',vgkq0l(:)
! reduced k+q vector: k'=k+q-K
  vkq0l(:)=vkq0l(:)-vgkq0l(:)
!-  write(*,*)'k+q0 in BZ, lattice coordinates:',vkq0l(:)
! search for index of reduced k+q vector 
  do jk=1,nkptnr
    if (r3taxi(vklnr(:,jk),vkq0l).lt.epslat) then
      k1(ik)=jk
      goto 10
    endif
  enddo
  write(*,*)
  write(*,'("Error(response): index of reduced k+q point is not found")')
  write(*,'(" Check q-vector coordinates")')
  write(*,*)
  stop
10 continue
enddo

! generate G+k vectors for entire BZ (this is required to compute wave-functions at each k-point)
do ik=1,nkptnr
  call gengpvec(vklnr(1,ik),vkcnr(1,ik),ngknr(ik),igkignr,vgklnr(1,1,ik),vgkcnr,gkcnr(1,ik),tpgkcnr(1,1,ik))
  call gensfacgp(ngknr(ik),vgkcnr,ngkmax,sfacgknr(1,1,ik))
enddo

! generate G+q0 vectors
do ig=1,ngvec
  vgq0c(:,ig)=vgc(:,ig)+vq0rc(:)
! get spherical coordinates and length of G+q0
  call sphcrd(vgq0c(:,ig),gq0(ig),tpgq0(:,ig))
! generate spherical harmonics for G+q0
  call genylm(lmaxvr,tpgq0(:,ig),ylmgq0(:,ig))
enddo

! generate structure factor for G+q0 vectors
call gensfacgp(ngvec,vgq0c,ngvec,sfacgq0)

! generate Bessel functions
do ig=1,ngvec
  do is=1,nspecies
    do ir=1,nrcmt(is)
!     |G+q0|*x
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
allocate(nnp(nkptnr,max_num_nnp,2))
allocate(docc(nkptnr,max_num_nnp))
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
        docc(ik,i1)=occsvnr(i,ik)-occsvnr(j,jk)
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

!- write(*,*)'size of wfmt arrays: ', &
!-  2*lmmaxvr*nrcmtmax*natmtot*nspinor*nstsv*16.d0/1024/1024,' Mb'
!- write(*,*)'size of wfir arrays: ', &
!-   2*ngrtot*nspinor*nstsv*16.d0/1024/1024,' Mb'

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

open(60,file='ZRHOFC.OUT',form='unformatted',status='replace')
write(60)nkptnr,ngvec_resp,max_num_nnp,igq0

write(50,*)
do ik=1,nkptnr
  write(*,'("k-point ",I4," out of ",I4)')ik,nkptnr
  call flush(50)
    
  jk = k1(ik)
  write(60)ik,jk
  write(60)num_nnp(ik)
  write(60)nnp(ik,1:num_nnp(ik),1:2)
  write(60)docc(ik,1:num_nnp(ik))

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
    call vnlrho(.true.,wfmt1(:,:,:,:,ist1),wfmt2(:,:,:,:,ist2),wfir1(:,:,ist1), &
      wfir2(:,:,ist2),zrhomt,zrhoir)
    call zrhoft(zrhomt,zrhoir,jlgq0r,ylmgq0,sfacgq0,ngvec_resp,zrhofc(1,i))
    write(60)zrhofc(1:ngvec_resp,i)
  enddo
enddo !ik
close(60)


deallocate(gshell,k1,vgq0c,gq0,tpgq0,sfacgq0,ylmgq0,jlgq0r,jl,vgklnr,   &
  vgkcnr,gkcnr,tpgkcnr,ngknr,sfacgknr,igkignr,occsvnr,num_nnp,nnp,docc, &
  zrhofc,evecfv,evecsv,apwalm,wfmt1,wfmt2,wfir1,wfir2,zrhomt,zrhoir)
endif !task.eq.400

if (task.eq.401) then
  write(50,'("Calculation of KS polarisability chi0")')
  nepts=400
  emax=1.d0
  
  open(60,file='ZRHOFC.OUT',form='unformatted',status='old')
  read(60)i1,ngvec_resp,max_num_nnp,igq0
  if (i1.ne.nkptnr) then
    write(*,*)'Error: k-mesh was changed'
    stop
  endif
  allocate(chi0(ngvec_resp,ngvec_resp,nepts))
  allocate(num_nnp(nkptnr))
  allocate(nnp(nkptnr,max_num_nnp,2))
  allocate(docc(nkptnr,max_num_nnp))
  allocate(evalsvnr(nstsv,nkptnr))
  allocate(w(nepts))
  allocate(zrhofc(ngvec_resp,max_num_nnp))  
  
  do i=1,nepts
    w(i)=dcmplx(1.d0*emax*i/nepts,0.02d0)
  enddo
  
  do ik=1,nkptnr
    call getevalsv(vklnr(1,ik),evalsvnr(1,ik))
  enddo
  
  chi0=dcmplx(0.d0,0.d0)
  do ik=1,nkptnr
    read(60)i1,jk
    if (i1.ne.ik) then
      write(*,*)'Error reading file ZRHOFC.OUT'
      stop
    endif
    read(60)num_nnp(ik)
    read(60)nnp(ik,1:num_nnp(ik),1:2)
    read(60)docc(ik,1:num_nnp(ik))
    
    do i=1,num_nnp(ik)
      read(60)zrhofc(1:ngvec_resp,i)
      do ie=1,nepts
        do ig=1,ngvec_resp
          do ig1=1,ngvec_resp
            chi0(ig,ig1,ie)=chi0(ig,ig1,ie)+docc(ik,i)/(evalsvnr(nnp(ik,i,1),ik)-evalsvnr(nnp(ik,i,2),jk)+w(ie))* &
              zrhofc(ig,i)*dconjg(zrhofc(ig1,i))
          enddo
        enddo
      enddo
    enddo
  enddo
  chi0=chi0/nkptnr
  close(60)
  open(60,file='resp.dat',form='formatted',status='replace')
  do ie=1,nepts
    write(60,*)ha2ev*dreal(w(ie)),dreal(chi0(igq0,igq0,ie)),dimag(chi0(igq0,igq0,ie))
  enddo
  close(60)
endif

close(50)
return
end
