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

real(8)                 :: vkq0l(3),t1
integer                 :: ik1,ik2,ist1,ist2,ik,jk,ig,ir,is,i,j,i1,i2,i3,ig1,ig2,ig3,ie
integer                 :: ivg1(3),ivg2(3)
complex(8), allocatable :: sfac3g(:)
real(8)                 :: jl(0:lmaxvr)
integer ,allocatable    :: ngknr(:)
integer ,allocatable    :: igkignr(:,:)
real(8) ,allocatable    :: vgklnr(:,:,:),vgkcnr(:,:),gkcnr(:,:),tpgkcnr(:,:,:)
complex(8) ,allocatable :: sfacgknr(:,:,:)
real(8) ,allocatable    :: occsvnr(:,:)
real                    :: cpu0,cpu1
complex(8) ,allocatable :: chi0(:,:,:)
complex(8) ,allocatable :: chi(:,:)
real(8) ,allocatable    :: evalsvnr(:,:)
complex(8) ,allocatable :: w(:)


! number of G-vectors for matrix elements calculation
integer                 :: ngvec_me
! q-vector in lattice coordinates
real(8)                 :: vq0l(3)
! q-vector in Cartesian coordinates
real(8)                 :: vq0c(3)
! reduced q-vector in lattice coordinates
real(8)                 :: vq0rl(3)
! reduced q-vector in Cartesian coordinates
real(8)                 :: vq0rc(3)
! G-vector which brings q to first BZ
integer                 :: vgq0l(3)
! index of G-vector which brings q to first BZ
integer                 :: igq0
! index of k'=k+q-K points
integer   , allocatable :: k1(:)
! indexes for fft-transform of u_{nk}^{*}u_{n'k'}exp{-iKx}, where k'=k+q-K
integer ,allocatable    :: igfft1(:,:)
! number of n,n' combinations of band indexes for each k-point
integer ,allocatable    :: num_nnp(:)
! maximum num_nnp over all k-points 
integer                 :: max_num_nnp
! pair of n,n' band indexes for each k-point
integer ,allocatable    :: nnp(:,:,:)
! array of Fourier coefficients of complex charge density
complex(8) ,allocatable :: zrhofc(:,:)
complex(8) ,allocatable :: zrhofc1(:,:)
! G+q vectors in Cart.coord.
real(8)   , allocatable :: vgq0c(:,:)
! length of G+q vectors
real(8)   , allocatable :: gq0(:)
! theta and phi angles of G+q vectors
real(8)   , allocatable :: tpgq0(:,:)
! sperical harmonics of G+q vectors
complex(8), allocatable :: ylmgq0(:,:)
! structure factor for G+q vectors
complex(8), allocatable :: sfacgq0(:,:)
! Bessel functions j_l(|G+q|x)
real(8) ,allocatable    :: jlgq0r(:,:,:,:)

! number of G-vectors for chi calculation
integer                 :: ngvec_chi

integer                 :: nkptnr_
integer                 :: ngsh_me_
real(8)                 :: vq0c_(3)

complex(8) ,allocatable :: mtrx1(:,:)

real(8) ,allocatable    :: vc(:)

integer                 :: nepts
real(8)                 :: emin,emax,de,eta

real(8) ,allocatable    :: docc(:,:)

real(8) ,allocatable    :: kkrel(:,:)

integer                 :: min_band_resp,max_band_resp
logical                 :: flg


real(8)                 :: norm1
complex(8)              :: znorm,zsum1(100)
complex(8)              :: znorm2(3)
real(8)                 :: v1(3)

integer                 :: info
integer ,allocatable    :: ipiv(:)
complex(8)              :: wt
real(8), parameter      :: au2ang=0.5291772108d0

! external functions
complex(8) ,external :: zfint
real(8) ,external       :: r3taxi


! initialise universal variables
call init0
call init1


if (task.eq.400) then
  allocate(vgklnr(3,ngkmax,nkptnr))
  allocate(vgkcnr(3,ngkmax))
  allocate(gkcnr(ngkmax,nkptnr))
  allocate(tpgkcnr(2,ngkmax,nkptnr))
  allocate(ngknr(nkptnr))
  allocate(sfacgknr(ngkmax,natmtot,nkptnr))
  allocate(igkignr(ngkmax,nkptnr))
! generate G+k vectors for entire BZ (this is required to compute 
!   wave-functions at each k-point)
  do ik=1,nkptnr
    call gengpvec(vklnr(1,ik),vkcnr(1,ik),ngknr(ik),igkignr(1,ik), &
      vgklnr(1,1,ik),vgkcnr,gkcnr(1,ik),tpgkcnr(1,1,ik))
    call gensfacgp(ngknr(ik),vgkcnr,ngkmax,sfacgknr(1,1,ik))
  enddo
  allocate(evecfv(nmatmax,nstfv))
  allocate(evecsv(nstsv,nstsv))
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
  allocate(wfmt1(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
  allocate(wfir1(ngrtot,nspinor,nstsv))
  allocate(wfmt2(lmmaxvr,nrcmtmax,natmtot,nspinor,nstsv))
  allocate(wfir2(ngrtot,nspinor,nstsv))
  allocate(zrhomt(lmmaxvr,nrcmtmax,natmtot))
  allocate(zrhoir(ngrtot))

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
endif

open(150,file='RESPONSE.OUT',form='formatted',status='replace')

if (task.eq.400) then
  write(150,'("Calculation of matrix elements &
    &<psi_{n,k}|e^{-i(G+q)x}|psi_{n'',k+q}>")')

! find number of G-vectors by given number of G-shells
  call getngvec(ngsh_me,ngvec_me)

  allocate(k1(nkptnr))
  allocate(vgq0c(3,ngvec))
  allocate(gq0(ngvec))
  allocate(tpgq0(2,ngvec))
  allocate(sfacgq0(ngvec,natmtot))
  allocate(sfac3g(natmtot))
  allocate(ylmgq0(lmmaxvr,ngvec)) 
  allocate(jlgq0r(nrcmtmax,0:lmaxvr,nspecies,ngvec))
  allocate(occsvnr(nstsv,nkptnr))
  allocate(igfft1(ngvec_me,nkptnr))

! q-vector in lattice coordinates
  do i=1,3
    vq0l(i)=1.d0*ivq0l(i)/ngridk(i)
  enddo

! find G-vector which brings q0 to first BZ
  vgq0l(:)=floor(vq0l(:))

! reduce q0 vector fo first BZ
  vq0rl(:)=vq0l(:)-vgq0l(:)

! check if we have enough G-shells to bring q-vector back to first BZ
  do ig=1,ngvec_me
    if (sum(abs(vgq0l(:)-ivg(:,ig))).eq.0) then
      igq0=ig
      goto 20
    endif
  enddo
  write(*,*)
  write(*,'("Error(response): not enough G-vectors to reduce q-vector &
    &to first BZ")')
  write(*,'(" Increase number of G-shells for matrix elements")')
  write(*,*)
  stop
20 continue

! get Cartesian coordinates of q-vector and reduced q-vector
  call r3mv(bvec,vq0l,vq0c)
  call r3mv(bvec,vq0rl,vq0rc)

  write(150,*)
  write(150,'("q-vector (lat.coord.)                        : ",&
    & 3G18.10)')vq0l
  write(150,'("q-vector (Cart.coord.) [a.u.]                : ",&
    & 3G18.10)')vq0c
  write(150,'("q-vector length [a.u.]                       : ",&
    & G18.10)')sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)
  write(150,'("q-vector length [1/A]                        : ",&
    & G18.10)')sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)/au2ang
  write(150,'("G-vector to reduce q to first BZ (lat.coord.): ",&
    & 3I4)')vgq0l
  write(150,'("index of G-vector                            : ",&
    & I4)')igq0
  write(150,'("reduced q-vector (lat.coord.)                : ",&
    & 3G18.10)')vq0rl
  write(150,'("reduced q-vector (Cart.coord.) [a.u.]        : ",&
    & 3G18.10)')vq0rc

! find k+q and reduce them to first BZ (this is required to utilize the 
!   periodical property of Bloch-states: |k>=|k+K>, where K is any vector 
!   of the reciprocal lattice)
  write(150,*)
  write(150,'(3X,"ik",10X,"k",19X,"k+q",16X,"K",12X,"k''=k+q-K",8X,"jk")')
  write(150,'(85("-"))')
  do ik=1,nkptnr
! k+q vector
    vkq0l(:)=vklnr(:,ik)+vq0rl(:)
! K vector
    ivg1(:)=floor(vkq0l(:))
! reduced k+q vector: k'=k+q-K
    vkq0l(:)=vkq0l(:)-ivg1(:)
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
10  continue
! search for new fft indexes
    do ig=1,ngvec_me
      ivg2(:)=ivg(:,ig)+ivg1(:)
      igfft1(ig,ik)=igfft(ivgig(ivg2(1),ivg2(2),ivg2(3)))
    enddo
    write(150,'(I4,2X,3F6.2,2X,3F6.2,2X,3I4,2X,3F6.2,2X,I4)') &
      ik,vklnr(:,ik),vkq0l+ivg1,ivg1,vkq0l,k1(ik)
  enddo

! get occupancy of states
  do ik=1,nkptnr
    call getoccsv(vklnr(1,ik),occsvnr(1,ik))
  enddo

  min_band_resp=1
  max_band_resp=nstsv
! setup n,n' stuff
! first, find the maximum size of nnp array
  max_num_nnp=0
  allocate(num_nnp(nkptnr))
  do ik=1,nkptnr
    jk=k1(ik)
    i1=0
    do i=min_band_resp,max_band_resp
      do j=min_band_resp,max_band_resp
        if (abs(occsvnr(i,ik)-occsvnr(j,jk)).gt.1d-10) i1=i1+1
      enddo
    enddo
    num_nnp(ik)=i1
    max_num_nnp=max(max_num_nnp,i1)
  enddo
  allocate(nnp(nkptnr,max_num_nnp,3))
  allocate(docc(nkptnr,max_num_nnp))
! second, setup the nnp array
  do ik=1,nkptnr
    jk=k1(ik)
    i1=0
    do i=min_band_resp,max_band_resp
      do j=min_band_resp,max_band_resp
        if (abs(occsvnr(i,ik)-occsvnr(j,jk)).gt.1d-10) then
          i1=i1+1
          nnp(ik,i1,1)=i
          nnp(ik,i1,2)=j
          docc(ik,i1)=occsvnr(i,ik)-occsvnr(j,jk)
        endif
      enddo
    enddo
  enddo

  allocate(zrhofc(ngvec_me,max_num_nnp))
  allocate(zrhofc1(ngvec_me,3))


  write(150,*)
! generate G+q' vectors, where q' is reduced q-vector
  do ig=1,ngvec_me
    vgq0c(:,ig)=vgc(:,ig)+vq0rc(:)
! get spherical coordinates and length of G+q'
    call sphcrd(vgq0c(:,ig),gq0(ig),tpgq0(:,ig))
! generate spherical harmonics for G+q'
    call genylm(lmaxvr,tpgq0(:,ig),ylmgq0(:,ig))
  enddo

! generate structure factor for G+q' vectors
  call gensfacgp(ngvec_me,vgq0c,ngvec,sfacgq0)
  
! generate Bessel functions j_l(|G+q'|x)
  do ig=1,ngvec_me
    do is=1,nspecies
      do ir=1,nrcmt(is)
        t1=gq0(ig)*rcmt(ir,is)
        call sbessel(lmaxvr,t1,jl)
        jlgq0r(ir,:,is,ig)=jl(:)
      enddo
    enddo
  enddo

  open(160,file='ZRHOFC.OUT',form='unformatted',status='replace')
  write(160)nkptnr,ngsh_me,ngvec_me,max_num_nnp,igq0
  write(160)vq0c(1:3)

  write(150,*)
  do ik=1,nkptnr
    write(150,'("k-point ",I4," out of ",I4)')ik,nkptnr
    
    jk=k1(ik)
  
    write(160)ik,jk
    write(160)num_nnp(ik)
    write(160)nnp(ik,1:num_nnp(ik),1:2)
    write(160)docc(ik,1:num_nnp(ik))
  
! generate wave-functions at k
    call getevecfv(vklnr(1,ik),vgklnr(:,:,ik),evecfv)
    call getevecsv(vklnr(1,ik),evecsv) 
    call match(ngknr(ik),gkcnr(:,ik),tpgkcnr(:,:,ik),sfacgknr(:,:,ik),apwalm)
    call genwfsv(.false.,ngknr(ik),igkignr(:,ik),evalsv(1,1),apwalm,evecfv, &
      evecsv,wfmt1,wfir1)

! test normalization    
    do i=1,nstsv
      call vnlrho(.true.,wfmt1(:,:,:,:,i),wfmt1(:,:,:,:,i),wfir1(:,:,i), &
        wfir1(:,:,i),zrhomt,zrhoir)
      znorm=zfint(zrhomt,zrhoir)
      if (abs(znorm-1.d0).gt.0.01d0) then
        write(150,'("Warning: bad norm ",G18.10," of wave-function ",&
          & I4," at k-point ",I4)')abs(znorm),i,ik
      endif
    enddo

! generate wave-functions at k'=k+q-K
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
      call zrhoft(zrhomt,zrhoir,jlgq0r,ylmgq0,sfacgq0,ngvec_me,igfft1(1,ik),zrhofc1)
      zrhofc(:,i)=zrhofc1(:,3)
    enddo
  
    write(160)zrhofc(1:ngvec_me,1:num_nnp(ik))
  enddo !ik
  close(160)


  deallocate(k1)
  deallocate(vgq0c)
  deallocate(gq0)
  deallocate(tpgq0)
  deallocate(sfacgq0)
  deallocate(ylmgq0)
  deallocate(jlgq0r)
  deallocate(vgklnr)
  deallocate(vgkcnr)
  deallocate(gkcnr)
  deallocate(tpgkcnr)
  deallocate(ngknr)
  deallocate(sfacgknr)
  deallocate(igkignr)
  deallocate(occsvnr)
  deallocate(num_nnp)
  deallocate(nnp)
  deallocate(docc)
  deallocate(zrhofc)
  deallocate(evecfv)
  deallocate(evecsv)
  deallocate(apwalm)
  deallocate(wfmt1)
  deallocate(wfmt2)
  deallocate(wfir1)
  deallocate(wfir2)
  deallocate(zrhomt)
  deallocate(zrhoir)

endif !task.eq.400

if (task.eq.401) then
  write(150,'("Calculation of KS polarisability chi0")')
  emin=0.d0
  emax=80.d0
  de=0.05d0
  eta=0.5d0
  
  write(150,*)'emin=',emin,' (eV)'
  write(150,*)'emax=',emax,' (eV)'
  write(150,*)'de=',de,' (eV)'
  write(150,*)'eta=',eta,' (eV)'
  
  
  
  nepts=1+(emax-emin)/de
  
  open(160,file='ZRHOFC.OUT',form='unformatted',status='old')
  read(160)nkptnr_,ngsh_me_,ngvec_me,max_num_nnp,igq0
  read(160)vq0c_(1:3)

  if (nkptnr_.ne.nkptnr) then
    write(*,*)
    write(*,'("Error(response): k-mesh was changed")')
    write(*,*)
    stop
  endif
  if (ngsh_me_.ne.ngsh_me) then
    write(*,*)
    write(*,'("Error(response): number of G-shells for matrix elements &
      &was changed")')
    write(*,*)
    stop
  endif
  if (ngsh_chi.gt.ngsh_me) then
    write(*,*)
    write(*,'("Error(response): wrong number of G-shells for calculation &
      &of chi")')
    write(*,*)
    stop
  endif

! find number of G-vectors by given number of G-shells
  call getngvec(ngsh_chi,ngvec_chi)
  
  if (igq0.gt.ngvec_chi) then
    write(*,*)
    write(*,'("Error(response): not enough G-vectors for calculation of &
      &chi")')
    write(*,'("  Increase number of G-shells for chi")')
    write(*,*)
    stop
  endif
    
  write(150,*)'igq0=',igq0
  
  allocate(chi0(ngvec_chi,ngvec_chi,nepts))
  allocate(chi(ngvec_chi,nepts))
  allocate(num_nnp(nkptnr))
  allocate(nnp(nkptnr,max_num_nnp,3))
  allocate(docc(nkptnr,max_num_nnp))
  allocate(evalsvnr(nstsv,nkptnr))
  allocate(w(nepts))
  allocate(zrhofc(ngvec_me,max_num_nnp))  
  
  allocate(vgq0c(3,ngvec_chi))
  allocate(vc(ngvec_chi))
  allocate(gq0(ngvec_chi))

  write(150,*)
  write(150,'("Coulomb potential matrix elements:")')
  write(150,'("ig    |G+q|    V")')
  write(150,'("----------------")')
  do ig=1,ngvec_chi
! generate G+q vectors  
    vgq0c(:,ig)=vgc(:,ig)+vq0c_(:)
    gq0(ig)=sqrt(vgq0c(1,ig)**2+vgq0c(2,ig)**2+vgq0c(3,ig)**2)
    vc(ig)=fourpi/gq0(ig)**2 
    write(150,'(I4,2x,2G18.10)')ig,gq0(ig),vc(ig)
  enddo
 
! setup energy mesh
  do i=1,nepts
    w(i)=dcmplx(emin+de*(i-1),eta)/ha2ev
  enddo

! get eigen-values  
  do ik=1,nkptnr
    call getevalsv(vklnr(1,ik),evalsvnr(1,ik))
  enddo
  
  chi0=dcmplx(0.d0,0.d0)
  do ik=1,nkptnr
    write(150,*)'ik=',ik,' out of ',nkptnr
    read(160)i1,jk
    if (i1.ne.ik) then
      write(*,*)'Error reading file ZRHOFC.OUT'
      stop
    endif
    read(160)num_nnp(ik)
    read(160)nnp(ik,1:num_nnp(ik),1:2)
    read(160)docc(ik,1:num_nnp(ik))
    read(160)zrhofc(1:ngvec_me,1:num_nnp(ik))
    
    do i=1,num_nnp(ik)
      do ie=1,nepts
        wt=docc(ik,i)/(evalsvnr(nnp(ik,i,1),ik)-evalsvnr(nnp(ik,i,2),jk)+w(ie))
        call zgerc(ngvec_chi,ngvec_chi,wt,zrhofc(1,i),1,zrhofc(1,i),1,chi0(1,1,ie),ngvec_chi)
      enddo !ie
    enddo !i
  enddo !ik
  close(160)
  
  chi0=chi0/nkptnr/omega
  
  open(160,file='chi0.dat',form='formatted',status='replace')
  do ie=1,nepts
    write(160,'(7G18.10)')dreal(w(ie)), &
      dreal(chi0(igq0,igq0,ie)),dimag(chi0(igq0,igq0,ie))
  enddo
  close(160)
  
  
  allocate(mtrx1(ngvec_chi0,ngvec_chi0))
    
  allocate(ipiv(ngvec_chi0))
  write(150,*)
  do ie=1,nepts
    write(150,*)'energy point ',ie,' out of ',nepts
! compute 1-chi0*V
    do i=1,ngvec_chi0
      do j=1,ngvec_chi0
        mtrx1(i,j)=-chi0(i,j,ie)*vc(j)
      enddo
      mtrx1(i,i)=dcmplx(1.d0,0.d0)+mtrx1(i,i)
    enddo
! solve [1-chi0*V]^{-1}*chi=chi0
    chi(:,ie)=chi0(:,igq0,ie)
    call zgesv(ngvec_chi,1,mtrx1,ngvec_chi,ipiv,chi(1,ie),ngvec_chi,info)
    if (info.ne.0) then
      write(*,*)'Error solving linear equations'
      write(*,*)'info=',info
      stop
    endif
  enddo !ie
  deallocate(ipiv)
  
  
  chi0=chi0/ha2ev/(au2ang)**3
  chi=chi/ha2ev/(au2ang)**3

!  allocate(kkrel(nepts,2))
!  kkrel=0.d0
!  
!  do ie=1,nepts
!    do i=1,nepts-1
!      if (i.eq.ie) cycle
!      kkrel(ie,1)=kkrel(ie,1)+2*(w(i+1)-w(i))*w(i)*dimag(chi0(igq0,igq0,i))/(w(i)**2-w(ie)**2)/pi
!      kkrel(ie,2)=kkrel(ie,2)-2*w(ie)*(w(i+1)-w(i))*dreal(chi0(igq0,igq0,i))/(w(i)**2-w(ie)**2)/pi
!    enddo
!  enddo
  
  open(160,file='resp.dat',form='formatted',status='replace')
  !do igq0=1,ngvec_chi0
    do ie=1,nepts
      write(160,'(7G18.10)')dreal(w(ie))*ha2ev, &
        dreal(chi0(igq0,igq0,ie)),dimag(chi0(igq0,igq0,ie)), &
        dreal(chi(igq0,ie)),dimag(chi(igq0,ie))
    enddo
    write(160,*)
  !enddo
  close(160)
endif

close(150)
return
end

subroutine getngvec(ngsh,ngv)
use modmain
implicit none
! arguments
integer ,intent(in)  :: ngsh
integer ,intent(out) :: ngv
! local variables
integer ,allocatable :: gshell(:)
integer              :: i

allocate(gshell(ngvec))

i=1
ngv=0
do while (i.le.ngsh)
  ngv=ngv+1
  gshell(ngv)=i
  if (abs(gc(ngv+1)-gc(ngv)).gt.epslat) then
    i=i+1
  endif
enddo 

write(150,*)
write(150,'("Number of G-shells  : ",I4)')ngsh
write(150,'("Number of G-vectors : ",I4)')ngv
write(150,*)
write(150,'("  G-vec.       lat.coord.      length(a.u.)   shell ")')
write(150,'(" ---------------------------------------------------")')
do i=1,ngv
  write(150,'(2X,I4,4X,3I5,4X,F12.6,5x,I4)')i,ivg(:,i),gc(i),gshell(i)
enddo

deallocate(gshell)

return
end