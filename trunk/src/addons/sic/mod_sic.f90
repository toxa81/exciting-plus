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
! number of poles on r-mesh
integer s_nrpole
! poles of r-mesh
real(8), allocatable :: s_rpole(:)
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
! spherical (theta,phi) coordinates of the covering points
real(8), allocatable :: s_tp(:,:)
! Cartesian coordinates of the covering points
real(8), allocatable :: s_x(:,:)
! weights for spherical integration
real(8), allocatable :: s_tpw(:)
! forward transformation from real spherical harmonics to coordinates
real(8), allocatable :: s_rlmf(:,:)
! backward transformation from coordinates to real spherical harmonics
real(8), allocatable :: s_rlmb(:,:)
! forward transformation from complex spherical harmonics to coordinates
complex(8), allocatable :: s_ylmf(:,:)
! backward transformation from coordinates to complex spherical harmonics
complex(8), allocatable :: s_ylmb(:,:)

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
complex(8), allocatable :: s_wvlm(:,:,:,:)

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
complex(8), external :: ylm_val
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
    zt1=exp(zi*dot_product(vkcnr(:,ik),vtc(:)))*wkptnr(ik)
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
    zt1=wkptnr(ik)*exp(zi*dot_product(x(:),vkcnr(:,ik)))/sqrt(omega)
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
integer ir1,lm,ir
real (8) x0,tp(2),dx
complex(8) zval 
complex(8) ylm(lmmaxwan)
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

subroutine s_spinor_func_val(x,flm,zval)
use modmain
implicit none
! arguments
real(8), intent(in) :: x(3)
complex(8), intent(in) :: flm(lmmaxwan,s_nr,nspinor)
complex(8), intent(out) :: zval(nspinor)
! local variables
integer ir1,lm,ir,ispn
real (8) x0,tp(2),dx
complex(8) zt1 
complex(8) ylm(lmmaxwan)
!
if (sum(x(:)**2).gt.(sic_wan_cutoff**2)) then
  zval=zzero
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
do ispn=1,nspinor
  zt1=zzero
  do lm=1,lmmaxwan
    zt1=zt1+(flm(lm,ir1,ispn)+dx*(flm(lm,ir1+1,ispn)-flm(lm,ir1,ispn)))*ylm(lm)
  enddo
  zval(ispn)=zt1
enddo
return
end subroutine

subroutine s_func_val2(x,f1lm,f2lm,zval1,zval2)
use modmain
implicit none
! arguments
real(8), intent(in) :: x(3)
complex(8), intent(in) :: f1lm(lmmaxwan,s_nr,nspinor)
complex(8), intent(in) :: f2lm(lmmaxwan,s_nr,nspinor)
complex(8), intent(out) :: zval1(nspinor)
complex(8), intent(out) :: zval2(nspinor)
! local variables
integer ir1,lm,ir,ispn
real (8) x0,tp(2),dx
complex(8) ylm(lmmaxwan)
complex(8) z1,z2
!
if (sum(x(:)**2).gt.(sic_wan_cutoff**2)) then
  zval1=zzero
  zval2=zzero
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
do ispn=1,nspinor
  z1=zzero
  z2=zzero
  do lm=1,lmmaxwan
    z1=z1+(f1lm(lm,ir1,ispn)+dx*(f1lm(lm,ir1+1,ispn)-f1lm(lm,ir1,ispn)))*ylm(lm)
    z2=z2+(f2lm(lm,ir1,ispn)+dx*(f2lm(lm,ir1+1,ispn)-f2lm(lm,ir1,ispn)))*ylm(lm)
  enddo
  zval1(ispn)=z1
  zval2(ispn)=z2
enddo
return
end subroutine 



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
complex(8), allocatable :: f1tp_(:,:)
complex(8), allocatable :: f2tp_(:,:)
complex(8), allocatable :: f2lm_(:,:)
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
  !allocate(f2lm_(lmmaxwan,s_nr))
  f2tp_=zzero
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(itp,x1,x2)
  do ir=1,s_nr
    do itp=1,s_ntp
      x1(:)=s_x(:,itp)*s_r(ir)
      x2(:)=pos1(:)+x1(:)-pos2(:)
      f2tp_(itp,ir)=s_func_val(x2,f2lm)
    enddo
  enddo !ir
!$OMP END PARALLEL DO
! convert to spherical coordinates
  call zgemm('T','N',s_ntp,s_nr,lmmaxwan,zone,s_ylmf,lmmaxwan,f1lm,&
    lmmaxwan,zzero,f1tp_,s_ntp)
  do ir=1,s_nr
    do itp=1,s_ntp
      zprod=zprod+dconjg(f1tp_(itp,ir))*f2tp_(itp,ir)*s_tpw(itp)*s_rw(ir)
    enddo
  enddo
! convert f2 to spherical harmonics
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

complex(8) function s_spinor_dot_ll(pos1,pos2,f1lm,f2lm)
use modmain
implicit none
! arguments
real(8), intent(in) :: pos1(3)
real(8), intent(in) :: pos2(3)
complex(8), intent(in) :: f1lm(lmmaxwan,s_nr,nspinor)
complex(8), intent(in) :: f2lm(lmmaxwan,s_nr,nspinor)
! local variables
complex(8), allocatable :: f1tp_(:,:)
complex(8), allocatable :: f2tp_(:,:,:)
complex(8), allocatable :: f2lm_(:,:)
complex(8) zprod
integer ir,itp,ispn
real(8) x1(3),x2(3)
complex(8), external :: zdotc
!
zprod=zzero
if (sum(abs(pos1-pos2)).lt.1d-10) then 
  do ispn=1,nspinor
    do ir=1,s_nr
      zprod=zprod+zdotc(lmmaxwan,f1lm(1,ir,ispn),1,f2lm(1,ir,ispn),1)*s_rw(ir)
    enddo
  enddo
else
  allocate(f1tp_(s_ntp,s_nr))
  allocate(f2tp_(s_ntp,s_nr,nspinor))
  !allocate(f2lm_(lmmaxwan,s_nr))
  f2tp_=zzero
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(itp,x1,x2)
  do ir=1,s_nr
    do itp=1,s_ntp
      x1(:)=s_x(:,itp)*s_r(ir)
      x2(:)=pos1(:)+x1(:)-pos2(:)
      call s_spinor_func_val(x2,f2lm,f2tp_(itp,ir,:))
    enddo
  enddo !ir
!$OMP END PARALLEL DO
  do ispn=1,nspinor
! convert to spherical coordinates
    call zgemm('T','N',s_ntp,s_nr,lmmaxwan,zone,s_ylmf,lmmaxwan,&
      f1lm(1,1,ispn),lmmaxwan,zzero,f1tp_,s_ntp)
    do ir=1,s_nr
      do itp=1,s_ntp
        zprod=zprod+dconjg(f1tp_(itp,ir))*f2tp_(itp,ir,ispn)*s_tpw(itp)*s_rw(ir)
      enddo
    enddo
  enddo
! convert f2 to spherical harmonics
!  call zgemm('T','N',lmmaxwan,s_nr,s_ntp,zone,s_ylmb,s_ntp,f2tp_,&
!    s_ntp,zzero,f2lm_,lmmaxwan)
!  do ir=1,s_nr
!    zprod=zprod+zdotc(lmmaxwan,f1lm(1,ir),1,f2lm_(1,ir),1)*s_rw(ir)
!  enddo
  deallocate(f1tp_,f2tp_)
endif
s_spinor_dot_ll=zprod
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
integer jr,ir,l,lm,ispn,lm1,lm2,lm3,lmmaxwanloc,lmloc
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
real(8), allocatable :: f2(:,:,:)
complex(8), allocatable :: f1(:,:,:)
complex(8), allocatable :: f3(:,:,:)
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
lmmaxwanloc=mpi_grid_map2(lmmaxwan,dims=(/dim_k,dim2/))
do lmloc=1,lmmaxwanloc
  lm=mpi_grid_map2(lmmaxwan,dims=(/dim_k,dim2/),loc=lmloc)
  l=lm2l(lm)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jr,t1)
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
!$OMP END PARALLEL DO
enddo
call mpi_grid_reduce(vhalm(1,1),lmmaxwan*s_nr,all=.true.)
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
  vxclm(:,:,ispn)=vxclm(:,:,ispn)+vhalm(:,:)
enddo
deallocate(rhotp,rholm,totrholm)
deallocate(vhalm,extp,ectp,exclm,vxtp,vctp)
! multiply Wannier function with potential and change sign
call timer_start(t_sic_wvprod)
allocate(f1(s_nr,lmmaxwan,nspinor))
allocate(f2(s_nr,lmmaxwan,nspinor))
allocate(f3(s_nr,lmmaxwan,nspinor))
! rearrange arrays in memory
do ispn=1,nspinor
  do lm=1,lmmaxwan
    f1(:,lm,ispn)=wanlm(lm,:,ispn)
    f2(:,lm,ispn)=vxclm(lm,:,ispn)
  enddo
enddo
f3=zzero
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(lm1,lm2,lm3,zt1,ispn)
do lmloc=1,lmmaxwanloc
  lm3=mpi_grid_map2(lmmaxwan,dims=(/dim_k,dim2/),loc=lmloc)
  do lm1=1,lmmaxwan
    do lm2=1,lmmaxwan
      zt1=gauntyry(lm2l(lm1),lm2l(lm2),lm2l(lm3),&
               lm2m(lm1),lm2m(lm2),lm2m(lm3))
      if (abs(zt1).gt.1d-12) then
        do ispn=1,nspinor
          !do ir=1,s_nr
          !  wvlm(lm3,ir,ispn)=wvlm(lm3,ir,ispn)-&
          !    wanlm(lm1,ir,ispn)*vxclm(lm2,ir,ispn)*zt1
          !enddo
          f3(:,lm3,ispn)=f3(:,lm3,ispn)-f1(:,lm1,ispn)*f2(:,lm2,ispn)*zt1
        enddo !ispn
      endif
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
call mpi_grid_reduce(f3(1,1,1),s_nr*lmmaxwan*nspinor,all=.true.) 
do ispn=1,nspinor
  do lm=1,lmmaxwan
    wvlm(lm,:,ispn)=f3(:,lm,ispn)
  enddo
enddo
deallocate(f1,f2,f3,vxclm)
call timer_stop(t_sic_wvprod)
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
