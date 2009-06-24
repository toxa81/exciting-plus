#ifdef _HDF5_
subroutine response_chi(ivq0m)
use modmain
implicit none
integer, intent(in) :: ivq0m(3)

integer igq0
! Kohn-Sham polarisability submatrix
complex(8), allocatable :: chi0s(:,:)
complex(8), allocatable :: chi0_GqGq(:)
! true polarisability
complex(8), allocatable :: chi(:)
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
complex(8), allocatable :: chi0_loc(:,:)
complex(8), allocatable :: lmbd(:,:)
real(8) fxca


integer ie,ig,i,j,ig1,ig2,ie1,ie2,idx0,bs
integer iv(3)
character*100 fname,qnm,path
logical, external :: root_cart

! we need Bcx and magnetization from STATE.OUT
if (lrtype.eq.1) call readstate

if (wproc) then
  write(150,*)
  if (lrtype.eq.0) then
    write(150,'("Calculation of charge polarizability chi")')  
  endif
  if (lrtype.eq.1) then
    write(150,'("Calculation of magnetic polarizability chi")')  
  endif
  call flushifc(150)
endif

call qname(ivq0m,qnm)
qnm="./"//trim(qnm)//"/"//trim(qnm)
fname=trim(qnm)//"_chi0.hdf5"
if (root_cart((/1,1,0/))) then
  call read_integer(nepts,1,trim(fname),'/parameters','nepts')
  call read_integer(lr_igq0,1,trim(fname),'/parameters','lr_igq0')
  call read_integer(gshme1,1,trim(fname),'/parameters','gshme1')
  call read_integer(gshme2,1,trim(fname),'/parameters','gshme2')
  call read_integer(gvecme1,1,trim(fname),'/parameters','gvecme1')
  call read_integer(gvecme2,1,trim(fname),'/parameters','gvecme2')
  call read_integer(ngvecme,1,trim(fname),'/parameters','ngvecme')
  call read_real8(vq0l,3,trim(fname),'/parameters','vq0l')
  call read_real8(vq0rl,3,trim(fname),'/parameters','vq0rl')
  call read_real8(vq0c,3,trim(fname),'/parameters','vq0c')
  call read_real8(vq0rc,3,trim(fname),'/parameters','vq0rc')
endif
call i_bcast_cart(comm_cart_110,nepts,1)
call i_bcast_cart(comm_cart_110,lr_igq0,1)
call i_bcast_cart(comm_cart_110,gshme1,1)
call i_bcast_cart(comm_cart_110,gshme2,1)
call i_bcast_cart(comm_cart_110,gvecme1,1)
call i_bcast_cart(comm_cart_110,gvecme2,1)
call i_bcast_cart(comm_cart_110,ngvecme,1)
call d_bcast_cart(comm_cart_110,vq0l,3)
call d_bcast_cart(comm_cart_110,vq0rl,3)
call d_bcast_cart(comm_cart_110,vq0c,3)
call d_bcast_cart(comm_cart_110,vq0rc,3)

if (wproc) then
  write(150,'("chi0 was calculated for ")')
  write(150,'("  G-shells  : ",I4," to ",I4)')gshme1,gshme2
  write(150,'("  G-vectors : ",I4," to ",I4)')gvecme1,gvecme2
  call flushifc(150)
endif

gshchi1=gshme1
gshchi2=gshme2
gvecchi1=gvecme1
gvecchi2=gvecme2

ngvecchi=gvecchi2-gvecchi1+1  
if (wproc) then 
  write(150,*)
  write(150,'("Minimum and maximum G-vectors for chi : ",2I4)')gvecchi1,gvecchi2
  write(150,'("Number of G-vectors : ",I4)')ngvecchi
  call flushifc(150)
endif

allocate(krnl(ngvecchi,ngvecchi))
allocate(ixcft(ngvec))
! construct kernel of the matrix equation
krnl=dcmplx(0.d0,0.d0)
! for charge response
if (lrtype.eq.0) then
  fxca=fxca0+fxca1*mpi_x(2)
  if (wproc) then
    write(150,*)
    write(150,'("fxc A : ",F8.4)')fxca
    write(150,*)
    write(150,'("Coulomb potential matrix elements:")')
    write(150,'("   ig        |G+q|        Vc-A/2  ")')
    write(150,'(" ------------------------------ ")')
  endif
  do ig=1,ngvecchi
! generate G+q vectors  
    vgq0c(:)=vgc(:,ig+gvecchi1-1)+vq0rc(:)
    gq0=sqrt(vgq0c(1)**2+vgq0c(2)**2+vgq0c(3)**2)
    krnl(ig,ig)=fourpi/gq0**2-fxca/2.d0
    if (wproc) then
      write(150,'(1X,I4,2X,2F12.6)')ig,gq0,abs(krnl(ig,ig))
    endif
  enddo
endif
! for magnetic response
if (lrtype.eq.1) then
  if (wproc) then
    write(150,*)
    write(150,'("Ixc matrix elements:")')
    write(150,'("   ig      |Ixc(G)|   ")')
    write(150,'(" -------------------- ")')
  endif
  call genixc(ixcft)
! contruct Ixc_{G,G'}=Ixc(G-G')
  do i=1,ngvecchi
    do j=1,ngvecchi
      ig1=gvecchi1+i-1
      ig2=gvecchi1+j-1
      iv(:)=-ivg(:,ig1)+ivg(:,ig2)
      krnl(i,j)=ixcft(ivgig(iv(1),iv(2),iv(3)))
    enddo
  enddo
  if (wproc) then
    do ig=1,2*ngvecchi
      write(150,'(1X,I4,2X,F12.6)')ig,abs(ixcft(ig))
    enddo
  endif
endif
if (wproc) then
  call flushifc(150)
endif


igq0=lr_igq0-gvecchi1+1

call idxbos(nepts,mpi_dims(1),mpi_x(1)+1,idx0,bs)
ie1=idx0+1
ie2=idx0+bs

allocate(lr_w(nepts))
allocate(chi0_loc(ngvecme,ngvecme))  
allocate(chi0s(ngvecchi,ngvecchi))
allocate(chi0_GqGq(nepts))
allocate(chi(nepts))
allocate(epsilon_GqGq(nepts))
allocate(epsilon_eff(nepts))
allocate(chi_scalar(nepts))
allocate(lmbd(ngvecchi,nepts))
lr_w=zzero
chi0_GqGq=zzero
chi=zzero
epsilon_GqGq=zzero
epsilon_eff=zzero
chi_scalar=zzero
lmbd=zzero

do ie=ie1,ie2
  if (root_cart((/0,1,0/))) then
#ifndef _PIO_
    do i=0,mpi_dims(1)-1
    do j=0,mpi_dims(3)-1
      if (mpi_x(1).eq.i.and.mpi_x(3).eq.j) then
#endif
        write(path,'("/iw/",I8.8)')ie
        call read_real8(lr_w(ie),2,trim(fname),trim(path),'w')
        call read_real8_array(chi0_loc,3,(/2,ngvecme,ngvecme/), &
          trim(fname),trim(path),'chi0')
#ifndef _PIO_      
      endif
      call barrier(comm_cart_101)
    enddo
    enddo
#endif
  endif
  call d_bcast_cart(comm_cart_010,lr_w(ie),2)
  call d_bcast_cart(comm_cart_010,chi0_loc,2*ngvecme*ngvecme)
! prepare chi0
  ig1=gvecchi1-gvecme1+1
  ig2=ig1+ngvecchi-1
  chi0s(1:ngvecchi,1:ngvecchi)=chi0_loc(ig1:ig2,ig1:ig2)
  chi0_GqGq(ie)=chi0s(igq0,igq0)
! solve matrix equation  
  call solve_chi(ngvecchi,igq0,chi0s,krnl,chi(ie), &
    chi_scalar(ie),epsilon_GqGq(ie),epsilon_eff(ie),lmbd(1,ie))
enddo

call d_reduce_cart(comm_cart_100,.false.,lr_w,2*nepts)
call d_reduce_cart(comm_cart_100,.false.,chi0_GqGq,2*nepts)
call d_reduce_cart(comm_cart_100,.false.,chi,2*nepts)
call d_reduce_cart(comm_cart_100,.false.,chi_scalar,2*nepts)
call d_reduce_cart(comm_cart_100,.false.,epsilon_GqGq,2*nepts)
call d_reduce_cart(comm_cart_100,.false.,epsilon_eff,2*nepts)
call d_reduce_cart(comm_cart_100,.false.,lmbd,2*ngvecchi*nepts)


if (root_cart((/1,0,0/))) then
  call write_chi(lr_igq0,ivq0m,chi0_GqGq,chi,chi_scalar, &
    epsilon_GqGq,epsilon_eff,lmbd)
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
      call dgemv('N',lmmaxvr,lmmaxvr,1.d0,rbshtapw,lmmaxapw,magmt(1,ir,ias,1),1, &
        0.d0,rftp2,1)
! calculate I(r)
      do itp=1,lmmaxvr
        if (abs(rftp2(itp)).lt.1d-10) rftp2(itp)=1d10
        zt1(itp)=dcmplx(rftp1(itp)/rftp2(itp),0.d0)
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
  if (abs(rt1).lt.1d-10) rt1=1d10
  zt3(ir)=dcmplx(bxcir(ir,1)/rt1,0.d0)*cfunir(ir)*omega
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


subroutine solve_chi(ngvecchi,igq0,chi0_in,krnl,chi,chi_scalar, &
  epsilon_GqGq,epsilon_eff,lmbd)
implicit none
integer, intent(in) :: ngvecchi
integer, intent(in) :: igq0
complex(8), intent(in) :: chi0_in(ngvecchi,ngvecchi)
complex(8), intent(in) :: krnl(ngvecchi,ngvecchi)
complex(8), intent(out) :: chi
complex(8), intent(out) :: chi_scalar
complex(8), intent(out) :: epsilon_GqGq
complex(8), intent(out) :: epsilon_eff
complex(8), intent(out) :: lmbd(ngvecchi)

complex(8), allocatable :: epsilon(:,:)
complex(8), allocatable :: chi_mtrx(:,:)
complex(8), allocatable :: mtrx1(:,:)
integer i

allocate(epsilon(ngvecchi,ngvecchi))
allocate(chi_mtrx(ngvecchi,ngvecchi))
allocate(mtrx1(ngvecchi,ngvecchi))
epsilon=dcmplx(0.d0,0.d0)
do i=1,ngvecchi
  epsilon(i,i)=dcmplx(1.d0,0.d0)
enddo
! compute epsilon=1-chi0*V
call zgemm('N','N',ngvecchi,ngvecchi,ngvecchi,dcmplx(-1.d0,0.d0), &
  chi0_in,ngvecchi,krnl,ngvecchi,dcmplx(1.d0,0.d0),epsilon, &
  ngvecchi)
mtrx1=epsilon
call diagzge(ngvecchi,mtrx1,lmbd) 
epsilon_GqGq=1.d0-chi0_in(igq0,igq0)*krnl(igq0,igq0)
chi_scalar=chi0_in(igq0,igq0)/epsilon_GqGq
call invzge(epsilon,ngvecchi)
epsilon_eff=1.d0/epsilon(igq0,igq0)
call zgemm('N','N',ngvecchi,ngvecchi,ngvecchi,dcmplx(1.d0,0.d0), &
  epsilon,ngvecchi,chi0_in,ngvecchi,dcmplx(0.d0,0.d0),chi_mtrx, &
  ngvecchi)
  chi=chi_mtrx(igq0,igq0)
deallocate(epsilon,chi_mtrx,mtrx1)
return
end

subroutine write_chi(igq0,ivq0m,chi0_in, &
  chi,chi_scalar,epsilon_GqGq,epsilon_eff,lmbd)
use modmain
implicit none
integer, intent(in) :: igq0
integer, intent(in) :: ivq0m(3)
complex(8), intent(in) :: chi0_in(nepts)
complex(8), intent(in) :: chi(nepts)
complex(8), intent(in) :: chi_scalar(nepts)
complex(8), intent(in) :: epsilon_GqGq(nepts)
complex(8), intent(in) :: epsilon_eff(nepts)
complex(8), intent(in) :: lmbd(ngvecchi,nepts)

real(8), allocatable :: func(:,:)
real(8) fxca
character*100 fname,qnm
character*10 c1,c2,c3,c4,c5
!character*4 c4
!character*2 c2
!character*1 c1
integer ie,i

call qname(ivq0m,qnm)
qnm="./"//trim(qnm)//"/"//trim(qnm)
fxca=fxca0+fxca1*mpi_x(2)
write(c2,'(F5.2)')fxca
write(c3,'(I8)')ngvecchi
write(c4,'(F6.3)')sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)/au2ang
write(c5,'(F5.3)')lr_eta

if (lrtype.eq.0) then
  fname=trim(qnm)//"__"//trim(adjustl(c4))//"__G_"//trim(adjustl(c3))//&
    "__A_"//trim(adjustl(c2))//"__.dat"
else
  fname=trim(qnm)//"_A"//c2//"_s"//c1//".dat"
endif
open(160,file=trim(fname),form='formatted',status='replace')
if (lrtype.eq.0) write(160,'("# charge density response")')
if (lrtype.eq.1) write(160,'("# magnetization density response")')
write(160,'("#")')
write(160,'("# Band interval (Ha) : ",2F8.2)')lr_e1,lr_e2
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
if (lrtype.eq.0) then
  write(160,'("#")')
  write(160,'("# fxc A : ",F8.4)')fxca
endif
write(160,'("#")')
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
close(160)

fname=trim(qnm)//"__"//trim(adjustl(c4))//"__G_"//trim(adjustl(c3))//&
  "__A_"//trim(adjustl(c2))//"_epseval__.dat"
open(160,file=trim(fname),form='formatted',status='replace')
do i=1,ngvecchi
  do ie=1,nepts
    write(160,'(4G18.10)')dreal(lr_w(ie))*ha2ev,abs(1.d0/lmbd(i,ie)),dreal(lmbd(i,ie)),dimag(lmbd(i,ie))  
  enddo
  write(160,'(" ")')
enddo
close(160)
return
end

#endif