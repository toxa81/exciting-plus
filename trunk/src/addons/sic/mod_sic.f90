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
!! vector -> index map
!  integer, allocatable :: ivtit(:,:,:)
!! translation limits along each lattice vector
!  integer tlim(2,3)
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
! coordinates of the unit vectors of the big sphere
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

! number of points for Lebedev mesh for LAPW muffin-tins 
integer mt_ntp
data mt_ntp/74/
! radial weights for LAPW muffin-tins
real(8), allocatable :: mt_rw(:,:)
! coordinates of the unit vectors of the MT-sphere
real(8), allocatable :: mt_spx(:,:)
! weights of Lebedev mesh for LAPW muffin-tins 
real(8), allocatable :: mt_tpw(:)
! forward transformation from complex spherical harmonics to coordinates
complex(8), allocatable :: mt_ylmf(:,:)


complex(8), allocatable :: s_wanlm(:,:,:,:)
!complex(8), allocatable :: s_pwanlm(:,:,:,:)
complex(8), allocatable :: s_wvlm(:,:,:,:)
!complex(8), allocatable :: s_pwvlm(:,:,:,:)

complex(8), allocatable :: s_wankmt(:,:,:,:,:,:)
complex(8), allocatable :: s_wankir(:,:,:,:)
complex(8), allocatable :: s_wvkmt(:,:,:,:,:,:)
complex(8), allocatable :: s_wvkir(:,:,:,:)

integer, parameter :: nwanprop=10
integer, parameter :: wp_normlm=1
integer, parameter :: wp_normtp=2
integer, parameter :: wp_rmswan=3
integer, parameter :: wp_rmsrho=4
integer, parameter :: wp_rmsrho13=5
integer, parameter :: wp_spread=6
integer, parameter :: wp_vha=7
integer, parameter :: wp_vxc=8
integer, parameter :: wp_vsic=9
integer, parameter :: wp_exc=10

contains

subroutine s_get_wffvval(x,ngp,vpc,vgpc,wffvmt,wffvit,wffvval)
use modmain
implicit none
real(8), intent(in) :: x(3)
integer, intent(in) :: ngp
real(8), intent(in) :: vpc(3)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(in) :: wffvmt(lmmaxvr,nufrmax,natmtot)
complex(8), intent(in) :: wffvit(nmatmax)
complex(8), intent(out) :: wffvval
integer is,ia,ias,ir0,io,l,j,i,lm,ig
integer ntr(3)
real(8) vtc(3),vr0(3),r0,tp(2),t1
real(8) ur(0:lmaxvr,nufrmax)
complex(8) zt1,ylm(lmmaxvr)
real(8) ya(nprad),c(nprad)
real(8), external :: polynom
logical, external :: vrinmt
!
wffvval=zzero
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
  do lm=1,lmmaxvr
    l=lm2l(lm)
    do io=1,nufr(l,is)
      wffvval=wffvval+wffvmt(lm,io,ias)*ur(l,io)*ylm(lm)
    enddo !io
  enddo !lm
  wffvval=wffvval*exp(zi*dot_product(vpc(:),vtc(:)))
else
  do ig=1,ngp
    zt1=exp(zi*dot_product(x(:),vgpc(:,ig)))/sqrt(omega)
    wffvval=wffvval+zt1*wffvit(ig)
  enddo
endif
return
end subroutine

!integer function s_gen_stepf(x)
!use modmain
!implicit none
!! arguments
!real(8), intent(in) :: x(3)
!! local variables
!integer is,ia,ir0
!integer ntr(3)
!real(8) vr0(3),r0
!logical, external :: vrinmt
!!
!s_gen_stepf=0
!if (.not.vrinmt(x,is,ia,ntr,vr0,ir0,r0)) s_gen_stepf=1
!return
!end function

subroutine s_get_ufrval(x,vpc,ias,ufrval)
use modmain
implicit none
! arguments
real(8), intent(in) :: x(3)
real(8), intent(in) :: vpc(3)
integer, intent(out) :: ias
complex(8), intent(out) :: ufrval(lmmaxvr,nufrmax)
! local variables
integer is,ia,ir0,io,l,j,i,lm
integer ntr(3)
real(8) vtc(3),vr0(3),r0,tp(2),t1
real(8) ur(0:lmaxvr,nufrmax)
complex(8) zt1,ylm(lmmaxvr)
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

subroutine s_get_wanval(twan,n,x,wanval)
use modmain
use mod_nrkp
use mod_wannier
implicit none
logical, intent(in) :: twan
integer, intent(in) :: n
real(8), intent(in) :: x(3)
complex(8), intent(out) :: wanval(nspinor)
!
integer is,ia,ias,ir0,io,l,j,i,lm,ig,ispn
integer ntr(3),ik,ikloc
real(8) x0(3),vtc(3),vr0(3),r0,tp(2),t1
real(8) ur(0:lmaxvr,nufrmax)
complex(8) zt1,zt2,ylm(lmmaxvr)
real(8) ya(nprad),c(nprad)
real(8), external :: polynom
logical, external :: vrinmt
complex(8) expigr(s_ngvec)
complex(8) zm1(lmmaxvr,nufrmax,nspinor),zm2(nspinor)
!
wanval=zzero
x0(:)=x(:)-wanpos(:,n)
if (sum(x0(:)**2).gt.(sic_wan_cutoff**2)) return
if (vrinmt(x,is,ia,ntr,vr0,ir0,r0).and.twan) then
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

!complex(8) function s_rinteg()
!use modmain
!implicit none
!complex(8), intent(in) :: fr(s_nr)
!real(8), allocatable :: fr1(:),fr2(:),gr(:),cf(:,:)
!
!allocate(fr1(s_nr))
!allocate(fr2(s_nr))
!allocate(gr(s_nr))
!allocate(cf(4,s_nr))
!zsummt=zzero
!do is=1,nspecies
!  do ia=1,natoms(is)
!    ias=idxas(ia,is)
!    do ir=1,nrmt(is)
!      zf1(ir)=zdotc(lmmaxvr,zfmt1(1,ir,ias),1,zfmt2(1,ir,ias),1)*spr(ir,is)**2
!      fr1(ir)=dreal(zf1(ir))
!      fr2(ir)=dimag(zf1(ir))
!    enddo
!    call fderiv(-1,nrmt(is),spr(1,is),fr1,gr,cf)
!    t1=gr(nrmt(is))
!    call fderiv(-1,nrmt(is),spr(1,is),fr2,gr,cf)
!    t2=gr(nrmt(is))
!    zsummt=zsummt+dcmplx(t1,t2)
!    !do ir=1,nrmt(is)-1
!    !  zsummt=zsummt+0.5d0*(zf1(ir)+zf1(ir+1))*(spr(ir+1,is)-spr(ir,is))
!    !enddo
!  end do
!end do
!

!return
!end

complex(8) function s_func_val(x,flm)
use modmain
implicit none
! arguments
real(8), intent(in) :: x(3)
complex(8), intent(in) :: flm(lmmaxwan,s_nr)
! local variables
integer ir1,lm,ir,np2,ir0,i,j
real (8) x0,tp(2),dx
complex(8) ylm(lmmaxwan),f3(lmmaxwan)
complex(8) zval 
real(8) c(nprad),f1(nprad),f2(nprad),t1,t2
real(8), external :: polynom
!
if (sum(x(:)**2).gt.(sic_wan_cutoff**2)) then
  s_func_val=zzero
  return
endif

call sphcrd(x,x0,tp)
call genylm(lmaxwan,tp,ylm)

!np2=nprad/2
!do ir=1,s_nr
!  if (s_r(ir).ge.x0) then
!    if (ir.le.np2) then
!      ir0=1
!    else if (ir.gt.(s_nr-np2)) then
!      ir0=s_nr-nprad+1
!    else
!      ir0=ir-np2
!    endif
!    x0=max(x0,s_r(1))
!    exit
!  endif
!enddo
!zval=zzero
!do lm=1,lmmaxwan
!  do j=1,nprad
!    i=ir0+j-1
!    f1(j)=dreal(flm(lm,i))
!    f2(j)=dimag(flm(lm,i))
!  enddo
!  t1=polynom(0,nprad,s_r(ir0),f1,c,x0)
!  t2=polynom(0,nprad,s_r(ir0),f2,c,x0)
!  zval=zval+dcmplx(t1,t2)*ylm(lm)
!enddo
!s_func_val=zval

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

subroutine s_func_plot1d(fname,np,p0,p1,p2,flm)
implicit none
character*(*), intent(in) :: fname
integer, intent(in) :: np
real(8), intent(in) :: p0(3)
real(8), intent(in) :: p1(3)
real(8), intent(in) :: p2(3)
complex(8), intent(in) :: flm(lmmaxwan,s_nr)

integer i
real(8) x(3),dx
complex(8) zt1
x(:)=(p2(:)-p1(:))/dble(np-1)
dx=sqrt(sum(x(:)**2))
open(220,file=trim(adjustl(fname)),form="formatted",status="replace")
do i=1,np
 x(:)=(p2(:)-p1(:))*(i-1)/dble(np-1)
 x(:)=x(:)-p0(:)
 zt1=s_func_val(x,flm)
 write(220,'(3G18.10)')dx*(i-1),dreal(zt1),dimag(zt1)
enddo
close(220)
return
end subroutine

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
!complex(8), allocatable :: f2lm_(:,:)
complex(8) zprod
integer ir,itp
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

subroutine s_gen_pot(wanlm,wantp,wvlm,wanprop)
use modmain
use modxcifc
implicit none
! arguments
complex(8), intent(in) :: wanlm(lmmaxwan,s_nr,nspinor)
complex(8), intent(in) :: wantp(s_ntp,s_nr,nspinor)
complex(8), intent(out) :: wvlm(lmmaxwan,s_nr,nspinor)
real(8), intent(out) :: wanprop(nwanprop)
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
! TODO: generalize for non-collinear case; vxc will become 2x2 matrix
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
wanprop(wp_vha)=0.d0
do ir=1,s_nr
  wanprop(wp_vha)=wanprop(wp_vha)+&
    ddot(lmmaxwan,totrholm(1,ir),1,vhalm(1,ir),1)*s_rw(ir)
enddo
! compute exc=<E_xc|rho>
wanprop(wp_exc)=0.d0
do ir=1,s_nr
  wanprop(wp_exc)=wanprop(wp_exc)+&
    ddot(lmmaxwan,totrholm(1,ir),1,exclm(1,ir),1)*s_rw(ir)
enddo
! compute vxc=<W_n|V_xc|W_n>; in the collinear case this is 
!  \sum_{\sigma} <V_xc^{\sigma}|rho${\sigma}>
wanprop(wp_vxc)=0.d0
do ispn=1,nspinor
  do ir=1,s_nr
    wanprop(wp_vxc)=wanprop(wp_vxc)+&
      ddot(lmmaxwan,rholm(1,ir,ispn),1,vxclm(1,ir,ispn),1)*s_rw(ir)
  enddo
enddo
! compute <V_n|rho>
wanprop(wp_vsic)=wanprop(wp_vha)+wanprop(wp_vxc)
! add Hartree potential to XC
do ispn=1,nspinor
  vxclm(:,:,ispn)=vhalm(:,:) !vxclm(:,:,ispn)+vhalm(:,:)
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

complex(8) function s_zfinp(tsh,tpw,ld,ng,zfmt1,zfmt2,zfir1,zfir2,zfrac)
use modmain
implicit none
logical, intent(in) :: tsh
logical, intent(in) :: tpw
integer, intent(in) :: ld
integer, intent(in) :: ng
complex(8), intent(in) :: zfmt1(ld,nrmtmax,natmtot)
complex(8), intent(in) :: zfmt2(ld,nrmtmax,natmtot)
complex(8), intent(in) :: zfir1(*)
complex(8), intent(in) :: zfir2(*)
complex(8), optional, intent(inout) :: zfrac(2)
!
integer is,ias,ir,itp,ig
complex(8) zsumir,zsummt,zt1
complex(8) zfr(nrmtmax) 
complex(8), external :: zdotc
! interstitial contribution
zsumir=zzero
if (tpw) then
  do ig=1,ng
    zsumir=zsumir+dconjg(zfir1(ig))*zfir2(ig)
  enddo
else
  do ir=1,ngrtot
    zsumir=zsumir+cfunir(ir)*conjg(zfir1(ir))*zfir2(ir)
  enddo
  zsumir=zsumir*omega/dble(ngrtot)
endif
! muffin-tin contribution
zsummt=zzero
do ias=1,natmtot
  is=ias2is(ias)
  do ir=1,nrmt(is)
    if (tsh) then
      zt1=zdotc(ld,zfmt1(1,ir,ias),1,zfmt2(1,ir,ias),1)
    else
      zt1=zzero
      do itp=1,mt_ntp
        zt1=zt1+mt_tpw(itp)*dconjg(zfmt1(itp,ir,ias))*zfmt2(itp,ir,ias)
      enddo
    endif
    zsummt=zsummt+zt1*mt_rw(ir,is)
  enddo
enddo !ias
s_zfinp=zsumir+zsummt
if (present(zfrac)) then
  zfrac(1)=zfrac(1)+zsummt
  zfrac(2)=zfrac(2)+zsumir
endif
return
end function























end module
