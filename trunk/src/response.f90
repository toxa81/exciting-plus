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

!- real(8)                 :: vq0l(3) = (/0.2d0,0.d0,0.d0/)
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
complex(8) ,allocatable :: chi(:,:,:)
real(8) ,allocatable    :: evalsvnr(:,:)
complex(8) ,allocatable :: w(:)
! number of G-shells for response
!- integer                 :: ngsh_resp
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

complex(8) ,allocatable :: mtrx1(:,:)

real(8) ,allocatable    :: vc(:,:)

integer                 :: nepts
real(8)                 :: emax

real(8) ,allocatable    :: docc(:,:)
integer    ,allocatable :: gshell(:)

real(8) ,allocatable    :: kkrel(:,:)

integer                 :: min_band_resp,max_band_resp
logical                 :: flg

real(8)                 :: norm1


! initialise universal variables
call init0
call init1

open(150,file='RESPONSE.OUT',form='formatted',status='replace')

if (task.eq.400) then
  write(150,'("Calculation of matrix elements <psi_{n,k}|e^{i(G+q)x}|psi_{n'',k+q}>")')

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

write(150,*)
write(150,'("Number of G-shells for response calculation  :",I4)')ngsh_resp
write(150,'("Number of G-vectors for response calculation :",I4)')ngvec_resp
write(150,*)
write(150,'("  G-vec.       lat.coord.      length(1/a.u.) shell")')
write(150,'(" ---------------------------------------------------")')
do ig=1,ngvec_resp
  write(150,'(2X,I4,4X,3I5,4X,F12.6,5x,I4)')ig,ivg(:,ig),gc(ig),gshell(ig)
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

write(150,*)
write(150,'("q-vector in lattice coordinates              : ",3G18.10)')vq0l
write(150,'("q-vector in Cartesian coordinates (1/a.u.)   : ",3G18.10)')vq0c
write(150,'("q-vector length (1/a.u.)                     : ",G18.10)')sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)
write(150,'("q-vector length (1/A)                        : ",G18.10)')sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)/0.529177d0
write(150,'("G-vector to reduce q to first BZ (lat.coord.): ",3I4)')vgq0l
write(150,'("Index of G-vector                            : ",I4)')igq0
write(150,'("Reduced q-vector (lat.coord.)                : ",3G18.10)')vq0rl

! find k+q and reduce them to first BZ (this is required to utilize the 
!   periodical property of Bloch-states: |k>=|k+K>, where K is any vector 
!   of the reciprocal lattice)
write(150,*)
write(150,'("  ik          k                   k+q                K            k''=k+q-K        jk")')
write(150,'("-------------------------------------------------------------------------------------")')
do ik=1,nkptnr
! k+q vector
  vkq0l(:)=vklnr(:,ik)+vq0l(:)
! K vector
  vgkq0l(:)=floor(vkq0l(:))
! reduced k+q vector: k'=k+q-K
  vkq0l(:)=vkq0l(:)-vgkq0l(:)
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
  write(150,'(I4,2X,3F6.2,2X,3F6.2,2X,3I4,2X,3F6.2,2X,I4)') &
    ik,vklnr(:,ik),vkq0l+vgkq0l,vgkq0l,vkq0l,k1(ik)
enddo

! generate G+k vectors for entire BZ (this is required to compute wave-functions at each k-point)
do ik=1,nkptnr
  call gengpvec(vklnr(1,ik),vkcnr(1,ik),ngknr(ik),igkignr(1,ik),vgklnr(1,1,ik),vgkcnr,gkcnr(1,ik),tpgkcnr(1,1,ik))
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
  !write(150,'(I4,2x,255F6.2)')ik,occsvnr(1:nstsv,ik)
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
      if (abs(occsvnr(i,ik)-occsvnr(j,jk)).gt.1d-8) i1=i1+1
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

open(160,file='ZRHOFC.OUT',form='unformatted',status='replace')
write(160)nkptnr,ngvec_resp,max_num_nnp,igq0
write(160)gq0(1:ngvec_resp)

write(150,*)
do ik=1,nkptnr
  write(150,'("k-point ",I4," out of ",I4)')ik,nkptnr
    
  jk = k1(ik)
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
    
! check normalization
  do i=1,nstsv
    call wfnorm(wfmt1(:,:,:,1,i),wfir1(:,1,i),norm1)
    if (abs(1.d0-norm1).gt.0.01) then
      write(*,*)
      write(*,'("Warning: norm of state ",I6," at k-point ",I6," is bad")')i,ik
      write(*,'(" norm: ",G18.10)')norm1
      write(*,'(" Try to turn off the symmetry or to decrease MT radii")')
      write(*,*)
    endif
    !call vnlrho(.true.,wfmt1(:,:,:,:,i),wfmt1(:,:,:,:,i),wfir1(:,:,i), &
    !  wfir1(:,:,i),zrhomt,zrhoir)
    ! for q=0 this should give the norm of wave-function
    !call zrhoft(zrhomt,zrhoir,jlgq0r,ylmgq0,sfacgq0,ngvec_resp,zrhofc(1,i))
    !write(150,*)i,zrhofc(1:ngvec_resp,i)
  enddo

! generate wave-functions at k+q
  call getevecfv(vklnr(1,jk),vgklnr(:,:,jk),evecfv)
  call getevecsv(vklnr(1,jk),evecsv) 
  call match(ngknr(jk),gkcnr(:,jk),tpgkcnr(:,:,jk),sfacgknr(:,:,jk),apwalm)
  call genwfsv(.false.,ngknr(jk),igkignr(:,jk),evalsv(1,1),apwalm,evecfv, &
    evecsv,wfmt2,wfir2)
  
  do i=1,num_nnp(ik)
    ist1=nnp(ik,i,1)
    ist2=nnp(ik,i,2)
    call vnlrho(.true.,wfmt1(:,:,:,:,ist2),wfmt2(:,:,:,:,ist1),wfir1(:,:,ist2), &
      wfir2(:,:,ist1),zrhomt,zrhoir)
    call zrhoft(zrhomt,zrhoir,jlgq0r,ylmgq0,sfacgq0,ngvec_resp,zrhofc(1,i))
        !write(150,'(4I3,255F12.6)')ik,jk,ist1,ist2,abs(zrhofc(1:ngvec_resp,i))
  enddo
  write(160)zrhofc(1:ngvec_resp,1:num_nnp(ik))

enddo !ik
close(160)


deallocate(gshell)
deallocate(k1)
deallocate(vgq0c)
deallocate(gq0)
deallocate(tpgq0)
deallocate(sfacgq0)
deallocate(ylmgq0)
deallocate(jlgq0r)
deallocate(jl)
deallocate(vgklnr)
deallocate(vgkcnr)
deallocate(gkcnr)
deallocate(tpgkcnr)
deallocate(ngknr)
stop
deallocate(sfacgknr,igkignr,occsvnr,num_nnp,nnp,docc, &
  zrhofc,evecfv,evecsv,apwalm,wfmt1,wfmt2,wfir1,wfir2,zrhomt,zrhoir)
endif !task.eq.400

if (task.eq.401) then
  write(150,'("Calculation of KS polarisability chi0")')
  nepts=2001
  emax=2.d0
  
  open(160,file='ZRHOFC.OUT',form='unformatted',status='old')
  read(160)i1,ngvec_resp,max_num_nnp,igq0
  if (i1.ne.nkptnr) then
    write(*,*)'Error: k-mesh was changed'
    stop
  endif
  allocate(chi0(ngvec_resp,ngvec_resp,nepts))
  allocate(chi(ngvec_resp,ngvec_resp,nepts))
  allocate(num_nnp(nkptnr))
  allocate(nnp(nkptnr,max_num_nnp,3))
  allocate(docc(nkptnr,max_num_nnp))
  allocate(evalsvnr(nstsv,nkptnr))
  allocate(w(nepts))
  allocate(zrhofc(ngvec_resp,max_num_nnp))  
  allocate(gq0(1:ngvec_resp))
  
  read(160)gq0(1:ngvec_resp)
  
  do i=1,nepts
    w(i)=dcmplx(emax*i/nepts,0.005d0)
  enddo
  
  do ik=1,nkptnr
    call getevalsv(vklnr(1,ik),evalsvnr(1,ik))
    write(150,'(I4,2x,255F6.2)')ik,evalsvnr(1:nstsv,ik)
  enddo
  
  chi0=dcmplx(0.d0,0.d0)
  do ik=1,nkptnr
    read(160)i1,jk
    if (i1.ne.ik) then
      write(*,*)'Error reading file ZRHOFC.OUT'
      stop
    endif
    read(160)num_nnp(ik)
    read(160)nnp(ik,1:num_nnp(ik),1:2)
    read(160)docc(ik,1:num_nnp(ik))
    read(160)zrhofc(1:ngvec_resp,1:num_nnp(ik))
    
    do i=1,num_nnp(ik)
      do ie=1,nepts
        do ig=1,ngvec_resp
          do ig1=1,ngvec_resp
            chi0(ig,ig1,ie)=chi0(ig,ig1,ie)+docc(ik,i)/(evalsvnr(nnp(ik,i,1),ik)-evalsvnr(nnp(ik,i,2),jk)+w(ie))* &
              zrhofc(ig,i)*dconjg(zrhofc(ig1,i))
          enddo
        enddo
      enddo
    enddo
  enddo !ik
  close(160)
  
  chi0=chi0/nkptnr/omega
  
  allocate(vc(ngvec_resp,ngvec_resp))
  allocate(mtrx1(ngvec_resp,ngvec_resp))
  vc=0.d0
  do ig=1,ngvec_resp
    vc(ig,ig)=fourpi/gq0(ig)**2
  enddo
  
  write(150,*)
  write(150,'("Coulomb potential matrix elements:")')
  do ig=1,ngvec_resp
    write(150,'(I4,2x,G18.10)')ig,vc(ig,ig)
  enddo
  
  
  do ie=1,nepts
! compute -chi0*V
    call zgemm('N','N',ngvec_resp,ngvec_resp,ngvec_resp,dcmplx(-1.d0,0.d0), &
      chi0(1,1,ie),ngvec_resp,vc,ngvec_resp,dcmplx(0.d0,0.d0),mtrx1,ngvec_resp)
! compute 1-chi0*V
    do i=1,ngvec_resp
      mtrx1(i,i)=mtrx1(i,i)+dcmplx(1.d0,0.d0)
    enddo
! compute [1-chi0*V]^{-1}
    call invzge(mtrx1,ngvec_resp)
! compute [1-chi0*V]^{-1}*chi0
    call zgemm('N','N',ngvec_resp,ngvec_resp,ngvec_resp,dcmplx(1.d0,0.d0), &
      mtrx1,ngvec_resp,chi0(1,1,ie),ngvec_resp,dcmplx(0.d0,0.d0),chi(1,1,ie),ngvec_resp)
  enddo
  
  
  chi0=(chi0/ha2ev)/(0.529)**3
  chi=(chi/ha2ev)/(0.529)**3
  
  allocate(kkrel(nepts,2))
  kkrel=0.d0
  
  do ie=1,nepts
    do i=1,nepts-1
      if (i.eq.ie) cycle
      kkrel(ie,1)=kkrel(ie,1)+2*(w(i+1)-w(i))*w(i)*dimag(chi0(igq0,igq0,i))/(w(i)**2-w(ie)**2)/pi
      kkrel(ie,2)=kkrel(ie,2)-2*w(ie)*(w(i+1)-w(i))*dreal(chi0(igq0,igq0,i))/(w(i)**2-w(ie)**2)/pi
    enddo
  enddo
  
  open(160,file='resp.dat',form='formatted',status='replace')
  do ie=1,nepts
    write(160,'(7G18.10)')dreal(w(ie)), &
      dreal(chi0(igq0,igq0,ie)),dimag(chi0(igq0,igq0,ie)), &
      kkrel(ie,1),kkrel(ie,2), &
      dreal(chi(igq0,igq0,ie)),dimag(chi(igq0,igq0,ie))
  enddo
  close(160)
endif

close(150)
return
end

subroutine wfnorm(wfmt,wfir,norm)
use modmain
implicit none
complex(8) ,intent(in)    :: wfmt(lmmaxvr,nrcmtmax,natmtot)
complex(8) ,intent(in)    :: wfir(ngrtot)
real(8)    ,intent(out  ) :: norm

integer is,nr,ia,ias,ir,l,m,lm
real(8) ,allocatable   :: fr1(:),gr(:),cf(:,:) 
real(8) t1,t2                                                                                           
complex(8) ,allocatable :: wfmt1(:,:)
real(8) sumir,summt,sum1
                                                                                                                                                  
allocate(fr1(nrcmtmax),gr(nrcmtmax),cf(3,nrcmtmax))
allocate(wfmt1(lmmaxvr,nrcmtmax))

! interstitial part
sumir=0.d0
do ir=1,ngrtot
  sumir=sumir+cfunir(ir)*abs(wfir(ir))**2
enddo
sumir=sumir*omega/ngrtot

! muffin-tin part
summt=0.d0
do is=1,nspecies
  nr=nrcmt(is)
  do ia=1,natoms(is) 
    ias=idxas(ia,is)
! transform to spherical harmonics
    call zgemm('N','N',lmmaxvr,nr,lmmaxvr,zone,zfshtvr,lmmaxvr,wfmt(1,1,ias), &
      lmmaxvr,zzero,wfmt1,lmmaxvr)
    do ir=1,nr
      sum1=0.d0
      do l=0,lmaxvr
        do m=-l,l
          lm=idxlm(l,m)
          sum1=sum1+abs(wfmt1(lm,ir))**2
	enddo !m
      enddo !l
      t1=rcmt(ir,is)**2
      fr1(ir)=sum1*t1
    enddo !ir
    call fderiv(-1,nr,rcmt(:,is),fr1,gr,cf)
    summt=summt+gr(nr)
  enddo !ia
enddo !is
  
norm=sumir+summt

deallocate(fr1,gr,cf)
deallocate(wfmt1)
return
end
