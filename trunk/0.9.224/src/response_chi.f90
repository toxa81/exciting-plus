subroutine response_chi(ivq0m)
use modmain
implicit none
integer, intent(in) :: ivq0m(3)

integer igq0
! Kohn-Sham polarisability submatrix
complex(8), allocatable :: chi0s(:,:,:,:)
complex(8), allocatable :: chi0_GqGq(:,:)
! true polarisability
complex(8), allocatable :: chi(:,:)
! G+q vector in Cartesian coordinates
real(8) vgq0c(3)
! length of G+q vector
real(8) gq0
! Fourier-transform of Ixc(r)=Bxc(r)/m(r)
complex(8), allocatable :: ixcft(:)

complex(8), allocatable :: epsilon_GqGq(:)
complex(8), allocatable :: chi_scalar(:)
complex(8), allocatable :: epsilon_eff(:)
! kernel of matrix equaiton
complex(8), allocatable :: krnl(:,:)


integer ie,ig,ngsh_me_,info,i,j,ig1,ig2,ispn
integer iv(3)
character*100 fname

! spin components : u(up) d(down)
! lrtype=0 : charge response
!   spin_me=1 : chi0_{uu}
!   spin_me=2 : chi0_{dd}
!   spin_me=3 : chi0_{uu} and chi0_{dd}
! lrtype=1 : magnetic response
!   spin_me=1 : chi0_{ud}
!   spin_me=2 : chi0_{du}
!   spin_me=3 : chi0_{ud} and chi0_{du}
!
! in the case spin_me=3 we can compute chi in four different ways:
! chi from chi0(1)
! chi from chi0(2)
! chi from chi0(1)+chi0(2)
! chi = chi(1)+chi(2)
!

! we need Bcx and magnetization from STATE.OUT
if (lrtype.eq.1) call readstate

if (iproc.eq.0) then

  write(150,*)
  if (lrtype.eq.0) then
    write(150,'("Calculation of charge polarizability chi")')  
  endif
  if (lrtype.eq.1) then
    write(150,'("Calculation of magnetic polarizability chi")')  
  endif
  if (lrtype.eq.2) then
    write(150,'("Calculation of chi = chi^{+-} + chi^{-+}")')  
  endif
  
  if (do_lr_io) then
    write(fname,'("CHI0[",I4.3,",",I4.3,",",I4.3,"].OUT")') &
      ivq0m(1),ivq0m(2),ivq0m(3)
    write(150,'("Reading file ",A40)')trim(fname)
    open(160,file=trim(fname),form='unformatted',status='old')
    read(160)nepts,lr_igq0
    read(160)gshme1,gshme2,gvecme1,gvecme2,ngvecme
    allocate(lr_w(nepts))
    read(160)lr_w(1:nepts)
    read(160)vq0l(1:3)
    read(160)vq0rl(1:3)
    read(160)vq0c(1:3)
    read(160)vq0rc(1:3)
    read(160)spin_me,nspin_chi0
    allocate(chi0(ngvecme,ngvecme,nepts,nspin_chi0))
    do ie=1,nepts
      read(160)chi0(1:ngvecme,1:ngvecme,ie,1:nspin_chi0)
    enddo
    if (task.eq.402) close(160)
    if (task.eq.403) close(160,status='delete')
  endif !do_lr_io
  write(150,'("chi0 was calculated for ")')
  write(150,'("  G-shells  : ",I4," to ",I4)')gshme1,gshme2
  write(150,'("  G-vectors : ",I4," to ",I4)')gvecme1,gvecme2
  if (spinpol.and.lrtype.eq.0) then 
    if (spin_me.eq.1) write(150,'("chi0 was calculated for spin up")')
    if (spin_me.eq.2) write(150,'("chi0 was calculated for spin dn")')
    if (spin_me.eq.3) write(150,'("chi0 was calculated for both spins")')
  endif
  if (spinpol.and.lrtype.eq.1) then 
    if (spin_me.eq.1) write(150,'("chi0 was calculated for up-dn")')
    if (spin_me.eq.2) write(150,'("chi0 was calculated for dn-up")')
    if (spin_me.eq.3) write(150,'("chi0 was calculated for up-dn and dn-up")')
  endif
  if (gvecchi1.lt.gvecme1.or.gvecchi1.gt.gvecme2) then
    write(150,*)
    write(150,'("Warning: minimum number of G-vectors was changed from ",&
      &I4," to ",I4)')gvecchi1,gvecme1
    gvecchi1=gvecme1 
  endif
  if (gvecchi2.lt.gvecme1.or.gvecchi2.gt.gvecme2) then
    write(150,*)
    write(150,'("Warning: maximum number of G-vectors was changed from ",&
      &I4," to ",I4)')gvecchi2,gvecme2
    gvecchi2=gvecme2 
  endif
  if (lr_igq0.lt.gvecchi1.or.lr_igq0.gt.gvecchi2) then
    write(*,*)
    write(*,'("Error(response_chi): not enough G-vectors for calculation of &
      &chi")')
    write(*,*)
    call pstop
  endif

  if (spinpol.and.afmchi0.and.nspin_chi0.eq.1.and.lrtype.eq.0) then
    write(150,'("AFM case: chi0 is multiplied by 2")')
    chi0=chi0*2.d0
  endif
  
  ngvecchi=gvecchi2-gvecchi1+1  
  write(150,*)
  write(150,'("Minimum and maximum G-vectors for chi : ",2I4)')gvecchi1,gvecchi2
  write(150,'("Number of G-vectors : ",I4)')ngvecchi

  igq0=lr_igq0-gvecchi1+1
  
  if (nspin_chi0.eq.1) then
    allocate(chi0s(ngvecchi,ngvecchi,nepts,1))
    allocate(chi0_GqGq(nepts,1))
    allocate(chi(nepts,1))
  else
    allocate(chi0s(ngvecchi,ngvecchi,nepts,3))
    allocate(chi0_GqGq(nepts,3))
    allocate(chi(nepts,4))
  endif
  allocate(epsilon_GqGq(nepts))
  allocate(epsilon_eff(nepts))
  allocate(chi_scalar(nepts))
  allocate(ixcft(ngvec))
  allocate(krnl(ngvecchi,ngvecchi))

! construct kernel of the matrix equation
  krnl=dcmplx(0.d0,0.d0)
! for charge response
  if (lrtype.eq.0) then
    write(150,*)
    write(150,'("Coulomb potential matrix elements:")')
    write(150,'("   ig        |G+q|        Vc    ")')
    write(150,'(" ------------------------------ ")')
    do ig=1,ngvecchi
! generate G+q vectors  
      vgq0c(:)=vgc(:,ig+gvecchi1-1)+vq0rc(:)
      gq0=sqrt(vgq0c(1)**2+vgq0c(2)**2+vgq0c(3)**2)
      krnl(ig,ig)=fourpi/gq0**2 
      write(150,'(1X,I4,2X,2F12.6)')ig,gq0,abs(krnl(ig,ig))
    enddo
  endif
! for magnetic response
  if (lrtype.eq.1) then
    call genixc(ixcft)
    write(150,*)
    write(150,'("Ixc matrix elements:")')
    write(150,'("   ig      |Ixc(G)|   ")')
    write(150,'(" -------------------- ")')
    do ig=1,2*ngvecchi
      write(150,'(1X,I4,2X,F12.6)')ig,abs(ixcft(ig))
    enddo
! contruct Ixc_{G,G'}=Ixc(G-G')
    do i=1,ngvecchi
      do j=1,ngvecchi
        ig1=gvecchi1+i-1
        ig2=gvecchi1+j-1
        iv(:)=-ivg(:,ig1)+ivg(:,ig2)
        krnl(i,j)=ixcft(ivgig(iv(1),iv(2),iv(3)))
      enddo
    enddo
  endif
! prepare chi0
  do ispn=1,nspin_chi0
    ig1=gvecchi1-gvecme1+1
    ig2=ig1+ngvecchi-1
    do ie=1,nepts
      chi0s(1:ngvecchi,1:ngvecchi,ie,ispn)=chi0(ig1:ig2,ig1:ig2,ie,ispn)
    enddo
  enddo
  if (nspin_chi0.eq.2) then
    chi0s(:,:,:,3)=chi0s(:,:,:,1)+chi0s(:,:,:,2)
  endif
  do ie=1,nepts
    chi0_GqGq(ie,:)=chi0s(igq0,igq0,ie,:)
  enddo
  
  if (nspin_chi0.eq.1) then
    call solve_chi(ngvecchi,nepts,igq0,chi0s(1,1,1,1),krnl,chi(1,1),chi_scalar, &
      epsilon_GqGq,epsilon_eff)
    call write_chi(lr_igq0,ivq0m,chi0_GqGq(1,1), &
      chi(1,1),chi_scalar,epsilon_GqGq,epsilon_eff,spin_me)
  else
! chi from chi0(1)
    call solve_chi(ngvecchi,nepts,igq0,chi0s(1,1,1,1),krnl,chi(1,1),chi_scalar, &
      epsilon_GqGq,epsilon_eff)
    call write_chi(lr_igq0,ivq0m,chi0_GqGq(1,1), &
      chi(1,1),chi_scalar,epsilon_GqGq,epsilon_eff,1)
! chi form chi0(2)
    call solve_chi(ngvecchi,nepts,igq0,chi0s(1,1,1,2),krnl,chi(1,2),chi_scalar, &
      epsilon_GqGq,epsilon_eff)
    call write_chi(lr_igq0,ivq0m,chi0_GqGq(1,2), &
      chi(1,2),chi_scalar,epsilon_GqGq,epsilon_eff,2)
! chi form chi0(1)+chi0(2)
    call solve_chi(ngvecchi,nepts,igq0,chi0s(1,1,1,3),krnl,chi(1,3),chi_scalar, &
      epsilon_GqGq,epsilon_eff)
    call write_chi(lr_igq0,ivq0m,chi0_GqGq(1,3), &
      chi(1,3),chi_scalar,epsilon_GqGq,epsilon_eff,3)
! chi = chi(1)+chi(2)    
    chi(:,4)=chi(:,1)+chi(:,2)
   call write_chi(lr_igq0,ivq0m,chi0_GqGq(1,3), &
      chi(1,4),chi_scalar,epsilon_GqGq,epsilon_eff,4)
  endif
      
  deallocate(lr_w)
  deallocate(chi0)
  deallocate(chi0s)
  deallocate(chi)
  deallocate(epsilon_GqGq)
  deallocate(epsilon_eff)
  deallocate(chi_scalar)
  deallocate(ixcft)
  deallocate(krnl)
  deallocate(chi0_GqGq)
  write(150,*)
  write(150,'("Done.")')

endif

return
end  

subroutine genixc(ixcft)
use modmain
implicit none
complex(8), intent(out) :: ixcft(ngvec)

integer is,ia,ias,l,m,lm,ir,ig,nr,itp
real(8) t1,t2,rt1
complex(8) zsum1,zsum2
real(8), allocatable :: rftp1(:)
real(8), allocatable :: rftp2(:)
real(8), allocatable :: jl(:,:)
real(8), allocatable :: fr1(:)
real(8), allocatable :: fr2(:)
real(8), allocatable :: gr(:)
real(8), allocatable :: cf(:,:)
complex(8), allocatable :: zt1(:)
complex(8), allocatable :: zt2(:,:,:)
complex(8), allocatable :: zt3(:)

allocate(rftp1(lmmaxvr))
allocate(rftp2(lmmaxvr))
allocate(jl(0:lmaxvr,nrmtmax))
allocate(fr1(nrmtmax))
allocate(fr2(nrmtmax))
allocate(gr(nrmtmax))
allocate(cf(3,nrmtmax))
allocate(zt1(lmmaxvr))
allocate(zt2(lmmaxvr,nrmtmax,natmtot))
allocate(zt3(ngrtot))

ixcft=dcmplx(0.d0,0.d0)

! calculate Ixc(r)=Bxc(r)/m(r) inside muffin-tins
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmt(is)
! transform Bxc and m from real spherical harmonics to spherical coordinates 
      call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtapw,lmmaxapw,bxcmt(1,ir,ias,1),1, &
        0.d0,rftp1,1)
!      call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtapw,lmmaxapw,vxcmt(1,ir,ias),1, &
!        0.d0,rftp1,1)
      call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtapw,lmmaxapw,magmt(1,ir,ias,1),1, &
        0.d0,rftp2,1)
! calculate I(r)
      do itp=1,lmmaxvr
        !if (abs(rftp2(itp)).lt.1d-10) rftp2(itp)=1d10
        zt1(itp)=dcmplx(rftp1(itp)/rftp2(itp),0.d0)
        !zt1(itp)=dcmplx(rftp2(itp),0.d0)
      enddo
! transform I(r) from spherical coordinates to complex spherical harmonics
      call zgemv('N',lmmaxvr,lmmaxvr,zone,zfshtvr,lmmaxvr,zt1,1,zzero, &
        zt2(1,ir,ias),1)
    enddo
  enddo
enddo
! calculate muffin-tin part of Fourier transform of Ixc(r) 
do ig=1,ngvec  
  do is=1,nspecies
    nr=nrmt(is)
! generate Bessel functions
    do ir=1,nr
      rt1=gc(ig)*spr(ir,is)
      call sbessel(lmaxvr,rt1,jl(0,ir))
    enddo
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ir=1,nr
        zsum1=dcmplx(0.d0,0.d0)
        do l=0,lmaxvr
          zsum2=dcmplx(0.d0,0.d0)
          do m=-l,l
            lm=idxlm(l,m)
            zsum2=zsum2+zt2(lm,ir,ias)*ylmg(lm,ig)
          enddo !m
          zsum1=zsum1+jl(l,ir)*dconjg(zil(l))*zsum2
        enddo !l
        rt1=spr(ir,is)**2
        fr1(ir)=dreal(zsum1)*rt1
        fr2(ir)=dimag(zsum1)*rt1
      enddo !ir
      call fderiv(-1,nr,spr(1,is),fr1,gr,cf)
      t1=gr(nr)
      call fderiv(-1,nr,spr(1,is),fr2,gr,cf)
      t2=gr(nr)
      ixcft(ig)=ixcft(ig)+fourpi*dconjg(sfacg(ig,ias))*dcmplx(t1,t2)
    enddo !ia
  enddo !is
enddo !ig


! calculate Ixc(r)=Bxc(r)/m(r) in interstitial
do ir=1,ngrtot
  rt1=magir(ir,1)
  !if (abs(rt1).lt.1d-10) rt1=1d10
  zt3(ir)=dcmplx(bxcir(ir,1)/rt1,0.d0)*cfunir(ir)*omega
!  zt3(ir)=dcmplx(vxcir(ir)/rt1,0.d0)*cfunir(ir)*omega
  !zt3(ir)=rt1*cfunir(ir)*omega
enddo       
call zfftifc(3,ngrid,-1,zt3)
do ig=1,ngvec
  ixcft(ig)=ixcft(ig)+zt3(igfft(ig))
enddo 
ixcft=ixcft/omega         

deallocate(rftp1)
deallocate(rftp2)
deallocate(jl)
deallocate(fr1)
deallocate(fr2)
deallocate(gr)
deallocate(cf)
deallocate(zt1)
deallocate(zt2)
deallocate(zt3)


return
end

subroutine solve_chi(ngvecchi,nepts,igq0,chi0_in,krnl,chi,chi_scalar, &
  epsilon_GqGq,epsilon_eff)
implicit none
integer, intent(in) :: ngvecchi
integer, intent(in) :: nepts
integer, intent(in) :: igq0
complex(8), intent(in) :: chi0_in(ngvecchi,ngvecchi,nepts)
complex(8), intent(in) :: krnl(ngvecchi,ngvecchi)
complex(8), intent(out) :: chi(nepts)
complex(8), intent(out) :: chi_scalar(nepts)
complex(8), intent(out) :: epsilon_GqGq(nepts)
complex(8), intent(out) :: epsilon_eff(nepts)

complex(8), allocatable :: epsilon(:,:)
complex(8), allocatable :: chi_mtrx(:,:)
integer i,ie

allocate(epsilon(ngvecchi,ngvecchi))
allocate(chi_mtrx(ngvecchi,ngvecchi))
do ie=1,nepts
  epsilon=dcmplx(0.d0,0.d0)
  do i=1,ngvecchi
    epsilon(i,i)=dcmplx(1.d0,0.d0)
  enddo
! compute epsilon=1-chi0*V
  call zgemm('N','N',ngvecchi,ngvecchi,ngvecchi,dcmplx(-1.d0,0.d0), &
    chi0_in(1,1,ie),ngvecchi,krnl,ngvecchi,dcmplx(1.d0,0.d0),epsilon, &
    ngvecchi)
!  epsilon_GqGq(ie)=epsilon(igq0,igq0)
  epsilon_GqGq(ie)=1.d0-chi0_in(igq0,igq0,ie)*krnl(igq0,igq0)
  chi_scalar(ie)=chi0_in(igq0,igq0,ie)/epsilon_GqGq(ie)
  call invzge(epsilon,ngvecchi)
  epsilon_eff(ie)=1.d0/epsilon(igq0,igq0)
  call zgemm('N','N',ngvecchi,ngvecchi,ngvecchi,dcmplx(1.d0,0.d0), &
    epsilon,ngvecchi,chi0_in(1,1,ie),ngvecchi,dcmplx(0.d0,0.d0),chi_mtrx, &
    ngvecchi)
  chi(ie)=chi_mtrx(igq0,igq0)
enddo !ie
deallocate(epsilon,chi_mtrx)
return
end

subroutine write_chi(igq0,ivq0m,chi0_in, &
  chi,chi_scalar,epsilon_GqGq,epsilon_eff,ispin_me)
use modmain
implicit none
integer, intent(in) :: igq0
integer, intent(in) :: ivq0m(3)
complex(8), intent(in) :: chi0_in(nepts)
complex(8), intent(in) :: chi(nepts)
complex(8), intent(in) :: chi_scalar(nepts)
complex(8), intent(in) :: epsilon_GqGq(nepts)
complex(8), intent(in) :: epsilon_eff(nepts)
integer, intent(in) :: ispin_me

real(8), allocatable :: func(:,:)
character*100 fname
integer ie

write(fname,'("response",I1,"[",I4.3,",",I4.3,",",I4.3,"].dat")') &
  ispin_me,ivq0m(1),ivq0m(2),ivq0m(3)
open(160,file=trim(fname),form='formatted',status='replace')
if (lrtype.eq.0) write(160,'("# charge density response")')
if (lrtype.eq.1) write(160,'("# magnetization density response")')

if (ispin_me.eq.1) write(160,'("# chi from chi0(1)")')
if (ispin_me.eq.2) write(160,'("# chi from chi0(2)")')
if (ispin_me.eq.3) write(160,'("# chi from chi0(1)+chi0(2)")')
if (ispin_me.eq.4) write(160,'("# chi = chi(1)+chi(2)")')
write(160,'("#")')
write(160,'("# k-mesh division                    : ",3I4)')ngridk(1),ngridk(2),ngridk(3)
write(160,'("# Energy mesh parameters             : ")')
write(160,'("#   maximum energy [eV]              : ", F9.4)')maxomega
write(160,'("#   energy step    [eV]              : ", F9.4)')domega
write(160,'("#   eta            [eV]              : ", F9.4)')lr_eta
write(160,'("# q-vector information               : ")')
write(160,'("#   q-vector (mesh coord.)           : ",3I4)')ivq0m
write(160,'("#   q-vector (lat. coord.)           : ",3F18.10)')vq0l
write(160,'("#   q-vector (Cart. coord.) [a.u.]   : ",3F18.10)')vq0c
write(160,'("#   q-vector length         [a.u.]   : ",3F18.10)')sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)
write(160,'("#   q-vector (Cart. coord.) [1/A]    : ",3F18.10)')vq0c/au2ang
write(160,'("#   q-vector length         [1/A]    : ",3F18.10)')sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)/au2ang
write(160,'("# G-vector information               : ")')
write(160,'("#   G-shells                         : ",2I4)')gshchi1,gshchi2
write(160,'("#   G-vectors                        : ",2I4)')gvecchi1,gvecchi2
write(160,'("#   index of Gq vector               : ",I4)')igq0
write(160,'("#")')
if (ispin_me.le.3) then
  write(160,'("# Definition of columns")')
  write(160,'("#   1: energy            [eV]")')
  write(160,'("#   2: -Re chi_0(Gq,Gq)  [1/eV/A^3]")')
  write(160,'("#   3: -Im chi_0(Gq,Gq)  [1/eV/A^3]")')
  write(160,'("#   4: -Re chi(Gq,Gq)    [1/eV/A^3]")')
  write(160,'("#   5: -Im chi(Gq,Gq)    [1/eV/A^3]")')
  write(160,'("#   6:  S(q,w)           [1/eV/A^3]")')
  write(160,'("#   7: -Re chi_scalar    [1/eV/A^3]")')
  write(160,'("#   8: -Im chi_scalar    [1/eV/A^3]")')
  write(160,'("#   9:  Re epsilon_eff             ")')
  write(160,'("#  10:  Im epsilon_eff             ")')
  write(160,'("#  11:  Re epsilon_GqGq            ")')
  write(160,'("#  12:  Im epsilon_GqGq            ")')
  write(160,'("#")')
  allocate(func(12,nepts))
  do ie=1,nepts
    func(1,ie)=dreal(lr_w(ie))*ha2ev
    func(2,ie)=-dreal(chi0_in(ie))/ha2ev/(au2ang)**3
    func(3,ie)=-dimag(chi0_in(ie))/ha2ev/(au2ang)**3
    func(4,ie)=-dreal(chi(ie))/ha2ev/(au2ang)**3
    func(5,ie)=-dimag(chi(ie))/ha2ev/(au2ang)**3
    func(6,ie)=-2.d0*dimag(chi(ie))/ha2ev/(au2ang)**3
    func(7,ie)=-dreal(chi_scalar(ie))/ha2ev/(au2ang)**3
    func(8,ie)=-dimag(chi_scalar(ie))/ha2ev/(au2ang)**3
    func(9,ie)=dreal(epsilon_eff(ie))
    func(10,ie)=dimag(epsilon_eff(ie))
    func(11,ie)=dreal(epsilon_GqGq(ie))
    func(12,ie)=dimag(epsilon_GqGq(ie))
    write(160,'(16F16.8)')func(1:12,ie)
  enddo
  deallocate(func)
else
  write(160,'("# Definition of columns")')
  write(160,'("#   1: energy            [eV]")')
  write(160,'("#   2: -Re chi(Gq,Gq)    [1/eV/A^3]")')
  write(160,'("#   3: -Im chi(Gq,Gq)    [1/eV/A^3]")')
  write(160,'("#   4:  S(q,w)           [1/eV/A^3]")')
  write(160,'("#")')
  allocate(func(4,nepts))
  do ie=1,nepts
    func(1,ie)=dreal(lr_w(ie))*ha2ev
    func(2,ie)=-dreal(chi(ie))/ha2ev/(au2ang)**3
    func(3,ie)=-dimag(chi(ie))/ha2ev/(au2ang)**3
    func(4,ie)=-2.d0*dimag(chi(ie))/ha2ev/(au2ang)**3
    write(160,'(16F16.8)')func(1:4,ie)
  enddo
  deallocate(func)
endif
close(160)
return
end

