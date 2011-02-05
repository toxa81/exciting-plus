module mod_sic
use mod_wannier

! matrix elements of Wannier potential <W_{n0}|V_{n0}|W_{n'T}>
complex(8), allocatable :: vwanme(:)
! LDA Hamiltonian in k-space in Wannier basis 
complex(8), allocatable :: sic_wann_h0k(:,:,:)
! LDA energies of Wannier functions
real(8), allocatable :: sic_wann_e0(:)
! number of SIC iterations
integer :: nsclsic
data nsclsic/3/
! current SIC iteration
integer :: isclsic
data isclsic/0/
! "DFT" energy before the SIC correction
real(8) :: engytot0
! total energy correction
real(8) :: sic_energy_tot
data sic_energy_tot/0.d0/
! potential contribution
real(8) :: sic_energy_pot
data sic_energy_pot/0.d0/
! kinetic contribution
real(8) :: sic_energy_kin
data sic_energy_kin/0.d0/
! cutoff distance for Wannier functions
real(8) :: sic_wan_cutoff
data sic_wan_cutoff/8.d0/
! cutoff distance for SIC marix elements <W_n|V_n|W_{n'T}>
real(8) :: sic_me_cutoff
data sic_me_cutoff/0.1d0/
! dot-product <W_{n\sigma}|f_{jk}> 
!  where f_{jk} is the first-variational Bloch state 
complex(8), allocatable :: sic_wb(:,:,:,:)
! dot-product <(W*V)_{n\sigma}|f_{jk}> 
!  where f_{jk} is the first-variational Bloch state 
complex(8), allocatable :: sic_wvb(:,:,:,:)
! dot-product <W_{n\sigma}|exp^{i(G+k)r}>
complex(8), allocatable :: sic_wgk(:,:,:,:)
! dot-product <(W*V)_{n\sigma}|exp^{i(G+k)r}>
complex(8), allocatable :: sic_wvgk(:,:,:,:)
! dot-product <W_{n\sigma}|u_{l}^{\alpha,\mu}*Y_{lm}>
complex(8), allocatable :: sic_wuy(:,:,:,:,:,:)
! dot-product <(W*V)_{n\sigma}|u_{l}^{\alpha,\mu}*Y_{lm}>
complex(8), allocatable :: sic_wvuy(:,:,:,:,:,:)

integer, allocatable :: sic_apply(:)
integer, allocatable :: sicw(:,:)

logical :: tsic_wv
data tsic_wv/.false./
logical :: tsic_arrays_allocated
data tsic_arrays_allocated/.false./

! maximum number of translation vectors
integer, parameter :: sic_maxvtl=1000

type(wannier_transitions) :: sic_wantran

type t_sic_orbitals
! total number of translations
  integer ntr
! translation vectors in lattice coordinates
  integer, allocatable :: vtl(:,:)
! translation vectors in Cartesian coordinates
  real(8), allocatable :: vtc(:,:)
! vector -> index map
  integer, allocatable :: ivtit(:,:,:)
! translation limits along each lattice vector
  integer tlim(2,3)
end type t_sic_orbitals

type(t_sic_orbitals) :: sic_orbitals

! maximum number of G-vectors for plane-wave expansion of Bloch functions
integer s_ngvec
! number of radial points in big spheres
integer s_nr
data s_nr/600/
! radial mesh of big spheres
real(8), allocatable :: s_r(:)
! weights for integration of radial functions
real(8), allocatable :: s_rw(:)
! maximum l for expansion in big spheres
integer lmaxwan
data lmaxwan/5/
! (lmaxwan+1)^2
integer lmmaxwan
! number of points on the big sphere
integer s_ntp
data s_ntp/266/
! coordinate of unit vectors on the sphere
real(8), allocatable :: s_spx(:,:)
! weights for spherical integration
real(8), allocatable :: s_tpw(:)
! forward transformation from real spherical harmonics to coordinates
real(8), allocatable :: s_rlmf(:,:)
! backward transformation from coordinates to real spherical harmonics
real(8), allocatable :: s_rlmb(:,:)
! backward transformation from coordinates to complex spherical harmonics
complex(8), allocatable :: s_ylmb(:,:)
! forward transformation from complex spherical harmonics to coordinates
complex(8), allocatable :: s_ylmf(:,:)


complex(8), allocatable :: s_wanlm(:,:,:,:)
complex(8), allocatable :: s_wvlm(:,:,:,:)
complex(8), allocatable :: s_wankmt(:,:,:,:,:,:)
complex(8), allocatable :: s_wankir(:,:,:,:)
complex(8), allocatable :: s_wvkmt(:,:,:,:,:,:)
complex(8), allocatable :: s_wvkir(:,:,:,:)

contains

!subroutine s_get_wffvval(ikloc,x,wffvmt,wffvit,wffvval)
!use modmain
!implicit none
!integer, intent(in) :: ikloc
!real(8), intent(in) :: x(3)
!complex(8), intent(in) :: wffvmt(lmmaxvr,nufrmax,natmtot)
!complex(8), intent(in) :: wffvit(nmatmax)
!complex(8), intent(out) :: wffvval
!integer is,ia,ias,ir0,io,l,j,i,lm,ig,ispn
!integer ntr(3),ik
!real(8) vrc0(3),vtc(3),vr0(3),r0,tp(2),t1
!real(8) ur(0:lmaxvr,nufrmax)
!complex(8) zt1,zt2,zt3,ylm(lmmaxvr)
!real(8) ya(nprad),c(nprad)
!real(8), external :: polynom
!logical, external :: vrinmt
!
!!ias1=wan_info(1,n1)
!!wanpos1(:)=atposc(:,ias2ia(ias1),ias2is(ias1))
!ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
!wffvval=zzero
!if (vrinmt(x,is,ia,ntr,vr0,ir0,r0)) then
!  ias=idxas(ia,is)
!  call sphcrd(vr0,t1,tp)
!  call genylm(lmaxvr,tp,ylm)
!  vtc(:)=ntr(1)*avec(:,1)+ntr(2)*avec(:,2)+ntr(3)*avec(:,3)
!  ur=0.d0
!  do l=0,lmaxvr
!    do io=1,nufr(l,is)
!      do j=1,nprad
!        i=ir0+j-1
!        ya(j)=ufr(i,l,io,ias2ic(ias))
!      end do
!      ur(l,io)=polynom(0,nprad,spr(ir0,is),ya,c,r0)
!    enddo !io
!  enddo !l
!  zt1=exp(zi*dot_product(vkc(:,ik),vtc(:)))
!  do lm=1,lmmaxvr
!    l=lm2l(lm)
!    do io=1,nufr(l,is)
!      wffvval=wffvval+wffvmt(lm,io,ias)*ur(l,io)*ylm(lm)
!    enddo !io
!  enddo !lm
!  wffvval=wffvval*zt1
!else
!  do ig=1,ngk(1,ik)
!    zt1=exp(zi*dot_product(x(:),vgkc(:,ig,1,ikloc)))/sqrt(omega)
!    wffvval=wffvval+zt1*wffvit(ig)
!  enddo
!endif
!return
!end subroutine


!subroutine s_get_pwval(ikloc,ig,x,pwval)
!use modmain
!implicit none
!integer, intent(in) :: ikloc
!integer, intent(in) :: ig
!real(8), intent(in) :: x(3)
!complex(8), intent(out) :: pwval
!
!integer is,ia,ir0
!integer ntr(3),ik
!real(8) vr0(3),r0
!logical, external :: vrinmt
!
!!ias1=wan_info(1,n1)
!!wanpos1(:)=atposc(:,ias2ia(ias1),ias2is(ias1))
!ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
!pwval=zzero
!if (.not.vrinmt(x,is,ia,ntr,vr0,ir0,r0)) then
!  pwval=exp(zi*dot_product(x(:),vgkc(:,ig,1,ikloc)))/sqrt(omega)
!endif
!return
!end subroutine

integer function s_gen_stepf(x)
use modmain
implicit none
! arguments
real(8), intent(in) :: x(3)
! local variables
integer is,ia,ir0
integer ntr(3),ik
real(8) vr0(3),r0
logical, external :: vrinmt
!
s_gen_stepf=0
if (.not.vrinmt(x,is,ia,ntr,vr0,ir0,r0)) s_gen_stepf=1
return
end function

!subroutine s_get_ufrval(ias,x,vpc,ufrval)
!use modmain
!implicit none
!! arguments
!integer, intent(in) :: ias
!real(8), intent(in) :: x(3)
!real(8), intent(in) :: vpc(3)
!complex(8), intent(out) :: ufrval(lmmaxvr,nufrmax)
!! local variables
!integer is,ia,ir0,io,l,j,i,lm,ig,ispn
!integer ntr(3)
!real(8) vrc0(3),vtc(3),vr0(3),r0,tp(2),t1
!real(8) ur(0:lmaxvr,nufrmax)
!complex(8) zt1,zt2,zt3,ylm(lmmaxvr)
!real(8) ya(nprad),c(nprad)
!real(8), external :: polynom
!logical, external :: vrinmt
!
!!ias1=wan_info(1,n1)
!!wanpos1(:)=atposc(:,ias2ia(ias1),ias2is(ias1))
!ufrval=zzero
!if (vrinmt(x,is,ia,ntr,vr0,ir0,r0)) then
!  if (idxas(ia,is).ne.ias) return
!  call sphcrd(vr0,t1,tp)
!  call genylm(lmaxvr,tp,ylm)
!  vtc(:)=ntr(1)*avec(:,1)+ntr(2)*avec(:,2)+ntr(3)*avec(:,3)
!  ur=0.d0
!  do l=0,lmaxvr
!    do io=1,nufr(l,is)
!      do j=1,nprad
!        i=ir0+j-1
!        ya(j)=ufr(i,l,io,ias2ic(ias))
!      end do
!      ur(l,io)=polynom(0,nprad,spr(ir0,is),ya,c,r0)
!    enddo !io
!  enddo !l
!  zt1=exp(zi*dot_product(vpc(:),vtc(:)))
!  do lm=1,lmmaxvr
!    l=lm2l(lm)
!    do io=1,nufr(l,is)
!      ufrval(lm,io)=ur(l,io)*ylm(lm)*zt1
!    enddo !io
!  enddo !lm
!endif
!return
!end subroutine

subroutine s_get_ufrval(x,vpc,ias,ufrval)
use modmain
implicit none
! arguments
real(8), intent(in) :: x(3)
real(8), intent(in) :: vpc(3)
integer, intent(out) :: ias
complex(8), intent(out) :: ufrval(lmmaxvr,nufrmax)
! local variables
integer is,ia,ir0,io,l,j,i,lm,ig,ispn
integer ntr(3)
real(8) vrc0(3),vtc(3),vr0(3),r0,tp(2),t1
real(8) ur(0:lmaxvr,nufrmax)
complex(8) zt1,zt2,zt3,ylm(lmmaxvr)
real(8) ya(nprad),c(nprad)
real(8), external :: polynom
logical, external :: vrinmt
!
ufrval=zzero
ias=-1
if (vrinmt(x,is,ia,ntr,vr0,ir0,r0)) then
  ias=idxas(ia,is)
  call sphcrd(vr0,t1,tp)
  call genylm(lmaxvr,tp,ylm)
  vtc(:)=ntr(1)*avec(:,1)+ntr(2)*avec(:,2)+ntr(3)*avec(:,3)
  ur=0.d0
  do l=0,lmaxvr
    do io=1,nufr(l,is)
      do j=1,nprad
        i=ir0+j-1
        ya(j)=ufr(i,l,io,ias2ic(ias))
      end do
      ur(l,io)=polynom(0,nprad,spr(ir0,is),ya,c,r0)
    enddo !io
  enddo !l
  zt1=exp(zi*dot_product(vpc(:),vtc(:)))
  do lm=1,lmmaxvr
    l=lm2l(lm)
    do io=1,nufr(l,is)
      ufrval(lm,io)=ur(l,io)*ylm(lm)*zt1
    enddo !io
  enddo !lm
endif
return
end subroutine




subroutine s_get_wanval(n,x,wanval)
use modmain
use mod_nrkp
use mod_wannier
implicit none
integer, intent(in) :: n
real(8), intent(in) :: x(3)
complex(8), intent(out) :: wanval(nspinor)
!
integer is,ia,ias,ir0,io,l,j,i,lm,ig,ispn
integer ntr(3),ik,ikloc
real(8) x0(3),vtc(3),vr0(3),r0,tp(2),t1
real(8) ur(0:lmaxvr,nufrmax)
complex(8) zt1,zt2,zt3,ylm(lmmaxvr)
real(8) ya(nprad),c(nprad)
real(8), external :: polynom
logical, external :: vrinmt
complex(8) expigr(s_ngvec)
complex(8) zm1(lmmaxvr,nufrmax,nspinor),zm2(nspinor)
!
wanval=zzero
x0(:)=x(:)-wanpos(:,n)
if (sum(x0(:)**2).gt.(sic_wan_cutoff**2)) return
if (vrinmt(x,is,ia,ntr,vr0,ir0,r0)) then
  ias=idxas(ia,is)
  call sphcrd(vr0,t1,tp)
  call genylm(lmaxvr,tp,ylm)
  vtc(:)=ntr(1)*avec(:,1)+ntr(2)*avec(:,2)+ntr(3)*avec(:,3)
  ur=0.d0
  do l=0,lmaxvr
    do io=1,nufr(l,is)
      do j=1,nprad
        i=ir0+j-1
        ya(j)=ufr(i,l,io,ias2ic(ias))
      end do
      ur(l,io)=polynom(0,nprad,spr(ir0,is),ya,c,r0)
    enddo !io
  enddo !l
  zm1=zzero
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    zt1=exp(zi*dot_product(vkcnr(:,ik),vtc(:)))/nkptnr
    do ispn=1,nspinor
      zm1(:,:,ispn)=zm1(:,:,ispn)+zt1*wann_unkmt(:,:,ias,ispn,n,ikloc)
    enddo
  enddo
  do ispn=1,nspinor
    do lm=1,lmmaxvr
      l=lm2l(lm)
      do io=1,nufr(l,is)
        wanval(ispn)=wanval(ispn)+zm1(lm,io,ispn)*ur(l,io)*ylm(lm)
      enddo !io
    enddo !lm
  enddo !ispn
else
  do ig=1,s_ngvec
    expigr(ig)=exp(zi*dot_product(x(:),vgc(:,ig)))  
  enddo
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    zt1=exp(zi*dot_product(x(:),vkcnr(:,ik)))/sqrt(omega)/nkptnr
    zm2=zzero
    do ig=1,ngknr(ikloc)
      zt2=expigr(igkignr(ig,ikloc))
      do ispn=1,nspinor
        zm2(ispn)=zm2(ispn)+zt2*wann_unkit(ig,ispn,n,ikloc)
      enddo
    enddo
    do ispn=1,nspinor
      wanval(ispn)=wanval(ispn)+zt1*zm2(ispn)
    enddo
  enddo
endif
return
end subroutine

complex(8) function s_func_val(x,flm)
use modmain
implicit none
! arguments
real(8), intent(in) :: x(3)
complex(8), intent(in) :: flm(lmmaxwan,s_nr)
! local variables
integer ir1,lm,itp,ir
real (8) x0,tp(2),dx
complex(8) ylm(lmmaxwan)
complex(8) zval 
!
if (sum(x(:)**2).gt.(sic_wan_cutoff**2)) then
  s_func_val=zzero
  return
endif

call sphcrd(x,x0,tp)
call genylm(lmaxwan,tp,ylm)
ir1=0
do ir=s_nr-1,1,-1
  if (s_r(ir).le.x0) then
    ir1=ir
    exit
  endif
enddo
if (ir1.eq.0) then
  ir1=1
  dx=0.d0
else
  dx=(x0-s_r(ir1))/(s_r(ir1+1)-s_r(ir1))
endif
zval=zzero
do lm=1,lmmaxwan
  zval=zval+(flm(lm,ir1)+dx*(flm(lm,ir1+1)-flm(lm,ir1)))*ylm(lm)
enddo
s_func_val=zval
return
end function

complex(8) function s_dot_ll(pos1,pos2,f1lm,f2lm)
use modmain
implicit none
! arguments
real(8), intent(in) :: pos1(3)
real(8), intent(in) :: pos2(3)
complex(8), intent(in) :: f1lm(lmmaxwan,s_nr)
complex(8), intent(in) :: f2lm(lmmaxwan,s_nr)
! local variables
complex(8), allocatable :: f2tp_(:,:)
complex(8), allocatable :: f1tp_(:,:)
complex(8), allocatable :: f2lm_(:,:)
complex(8) zprod
integer ias1,ias2,ir,itp,lm
real(8) x1(3),x2(3)
complex(8), external :: zdotc
!
zprod=zzero
if (sum(abs(pos1-pos2)).lt.1d-10) then 
  do ir=1,s_nr
    zprod=zprod+zdotc(lmmaxwan,f1lm(1,ir),1,f2lm(1,ir),1)*s_rw(ir)
  enddo
else
  allocate(f1tp_(s_ntp,s_nr))
  allocate(f2tp_(s_ntp,s_nr))
!  allocate(f2lm_(lmmaxwan,s_nr))
  f2tp_=zzero
  do ir=1,s_nr
    do itp=1,s_ntp
      x1(:)=s_spx(:,itp)*s_r(ir)
      x2(:)=pos1(:)+x1(:)-pos2(:)
      f2tp_(itp,ir)=s_func_val(x2,f2lm)
    enddo
  enddo !ir
  call zgemm('T','N',s_ntp,s_nr,lmmaxwan,zone,s_ylmf,lmmaxwan,f1lm,&
    lmmaxwan,zzero,f1tp_,s_ntp)
  do ir=1,s_nr
    do itp=1,s_ntp
      zprod=zprod+dconjg(f1tp_(itp,ir))*f2tp_(itp,ir)*s_tpw(itp)*s_rw(ir)
    enddo
  enddo
!  call zgemm('T','N',lmmaxwan,s_nr,s_ntp,zone,s_ylmb,s_ntp,f2tp_,&
!    s_ntp,zzero,f2lm_,lmmaxwan)
!  do ir=1,s_nr
!    zprod=zprod+zdotc(lmmaxwan,f1lm(1,ir),1,f2lm_(1,ir),1)*s_rw(ir)
!  enddo
  deallocate(f1tp_,f2tp_)
endif
s_dot_ll=zprod
return
end function

subroutine s_gen_pot(wanlm,wantp,wvlm,vha,vxc,vsic,exc)
use modmain
use modxcifc
implicit none
! arguments
complex(8), intent(in) :: wanlm(lmmaxwan,s_nr,nspinor)
complex(8), intent(in) :: wantp(s_ntp,s_nr,nspinor)
complex(8), intent(out) :: wvlm(lmmaxwan,s_nr,nspinor)
real(8), intent(out) :: vha
real(8), intent(out) :: vxc
real(8), intent(out) :: vsic
real(8), intent(out) :: exc
! local variables
integer jr,ir,l,lm,ispn,lm1,lm2,lm3
real(8) t1
complex(8) zt1
real(8), allocatable :: rhotp(:,:,:)
real(8), allocatable :: rholm(:,:,:)
real(8), allocatable :: totrholm(:,:)
real(8), allocatable :: vhalm(:,:)
real(8), allocatable :: extp(:,:)
real(8), allocatable :: ectp(:,:)
real(8), allocatable :: vxtp(:,:,:)
real(8), allocatable :: vctp(:,:,:)
real(8), allocatable :: exclm(:,:)
real(8), allocatable :: vxclm(:,:,:)
real(8), external :: ddot
complex(8), external :: gauntyry
!
allocate(rhotp(s_ntp,s_nr,nspinor))
allocate(rholm(lmmaxwan,s_nr,nspinor))
allocate(totrholm(lmmaxwan,s_nr))
allocate(vhalm(lmmaxwan,s_nr))
allocate(extp(s_ntp,s_nr))
allocate(ectp(s_ntp,s_nr))
allocate(exclm(lmmaxwan,s_nr))
allocate(vxtp(s_ntp,s_nr,nspinor))
allocate(vctp(s_ntp,s_nr,nspinor))
allocate(vxclm(lmmaxwan,s_nr,nspinor))

totrholm=0.d0
do ispn=1,nspinor
  rhotp(:,:,ispn)=abs(wantp(:,:,ispn))**2
! convert spin density to real spherical harmonic expansion
  call dgemm('T','N',lmmaxwan,s_nr,s_ntp,1.d0,s_rlmb,s_ntp,rhotp(1,1,ispn),&
    s_ntp,0.d0,rholm(1,1,ispn),lmmaxwan)
  totrholm(:,:)=totrholm(:,:)+rholm(:,:,ispn)
enddo
! compute Hartree potential
vhalm=0.d0
do lm=1,lmmaxwan
  l=lm2l(lm)
  do ir=1,s_nr
    t1=0.d0
    do jr=1,ir
      t1=t1+(s_r(jr)**l/s_r(ir)**(l+1))*totrholm(lm,jr)*s_rw(jr)
    enddo
    do jr=ir+1,s_nr
      t1=t1+(s_r(ir)**l/s_r(jr)**(l+1))*totrholm(lm,jr)*s_rw(jr)
    enddo
    vhalm(lm,ir)=t1*fourpi/(2*l+1)
  enddo
enddo
! compute XC potential
if (spinpol) then
  call xcifc(xctype,n=s_ntp*s_nr,rhoup=rhotp(1,1,1),rhodn=rhotp(1,1,2),&
    ex=extp,ec=ectp,vxup=vxtp(1,1,1),vxdn=vxtp(1,1,2),vcup=vctp(1,1,1),&
    vcdn=vctp(1,1,2))
else
  call xcifc(xctype,n=s_ntp*s_nr,rho=rhotp,ex=extp,ec=ectp,vx=vxtp,vc=vctp)
endif
! save XC energy density in extp
extp(:,:)=extp(:,:)+ectp(:,:)
! expand in real spherical harmonics
call dgemm('T','N',lmmaxwan,s_nr,s_ntp,1.d0,s_rlmb,s_ntp,extp,s_ntp,0.d0,&
  exclm,lmmaxwan)
! save XC potential in vxtp and expand in real spherical harmonics   
do ispn=1,nspinor
  vxtp(:,:,ispn)=vxtp(:,:,ispn)+vctp(:,:,ispn)
  call dgemm('T','N',lmmaxwan,s_nr,s_ntp,1.d0,s_rlmb,s_ntp,vxtp(1,1,ispn),&
    s_ntp,0.d0,vxclm(1,1,ispn),lmmaxwan)
enddo
! compute vha=<V_h|rho>
vha=0.d0
do ir=1,s_nr
  vha=vha+ddot(lmmaxwan,totrholm(1,ir),1,vhalm(1,ir),1)*s_rw(ir)
enddo
! compute exc=<E_xc|rho>
exc=0.d0
do ir=1,s_nr
  exc=exc+ddot(lmmaxwan,totrholm(1,ir),1,exclm(1,ir),1)*s_rw(ir)
enddo
! compute vxc=<V_xc|rho>
vxc=0.d0
do ispn=1,nspinor
  do ir=1,s_nr
    vxc=vxc+ddot(lmmaxwan,rholm(1,ir,ispn),1,vxclm(1,ir,ispn),1)*s_rw(ir)
  enddo
enddo
! compute <V_n|rho>
vsic=vha+vxc
! add Hartree potential to XC
do ispn=1,nspinor
  vxclm(:,:,ispn)=vxclm(:,:,ispn)+vhalm(:,:)
enddo
! multiply Wannier function with potential and change sign
wvlm=zzero
do lm1=1,lmmaxwan
  do lm2=1,lmmaxwan
    do lm3=1,lmmaxwan
      zt1=gauntyry(lm2l(lm1),lm2l(lm2),lm2l(lm3),&
               lm2m(lm1),lm2m(lm2),lm2m(lm3))
      if (abs(zt1).gt.1d-12) then
        do ispn=1,nspinor
          do ir=1,s_nr
            wvlm(lm3,ir,ispn)=wvlm(lm3,ir,ispn)-&
              wanlm(lm1,ir,ispn)*vxclm(lm2,ir,ispn)*zt1
          enddo
        enddo !ispn
      endif
    enddo
  enddo
enddo
deallocate(rhotp,rholm,totrholm)
deallocate(vhalm,extp,ectp,exclm,vxtp,vctp,vxclm)
end subroutine     

end module
