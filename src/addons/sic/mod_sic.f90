module mod_sic
use mod_wannier
!---------------------!
! spherical expansion !
!---------------------!
! cutoff distance for Wannier functions (probaly will be removed)
!real(8) sic_wan_cutoff
!data sic_wan_cutoff/1.d6/



! radius of the big sphere
real(8) s_rmax
data s_rmax/1.d10/
! innter ("safe") radius used to compute radial integrals 
real(8) s_rmin
! outer radius of the Wigner-Seitz cell of the macrocrystal
real(8) sic_wan_rwsmax
! inner radius of the Wigner-Seitz cell of the macrocrystal
real(8) sic_wan_rwsmin
! number of radial points
integer s_nr
data s_nr/600/
! number of radial points for sic_wan_rmin
integer s_nr_min
! radial mesh of big spheres
real(8), allocatable :: s_r(:)
! number of poles of the radial mesh
integer s_nrpole
! poles of r-mesh
real(8), allocatable :: s_rpole(:)
! radial integration weights
real(8), allocatable :: s_rw(:)
! maximum l for expansion in big spheres
integer lmaxwan
data lmaxwan/10/
! (lmaxwan+1)^2
integer lmmaxwan
! number of covering spherical points
integer s_ntp
! spherical (theta,phi) coordinates of the covering points
real(8), allocatable :: s_tp(:,:)
! Cartesian coordinates of the covering points
real(8), allocatable :: s_x(:,:)
! weights for spherical integration
real(8), allocatable :: s_tpw(:)
! TODO: "forward" <=> "backward" convention
! forward transformation from real spherical harmonics to coordinates
real(8), allocatable :: s_rlmf(:,:)
! backward transformation from coordinates to real spherical harmonics
real(8), allocatable :: s_rlmb(:,:)
! forward transformation from complex spherical harmonics to coordinates
complex(8), allocatable :: s_ylmf(:,:)
! backward transformation from coordinates to complex spherical harmonics
complex(8), allocatable :: s_ylmb(:,:)
! number of iterations for recursive reconstruction of expansion coefficients
integer sic_bsht_niter
data sic_bsht_niter/10/
! controls number of spherical points for different types of spherical meshes 
integer sic_smesh_n
data sic_smesh_n/11/
! Ylm expansion coefficients of Wannier functions
complex(8), allocatable :: s_wlm(:,:,:,:)
! Ylm expansion coefficients of W*V, W is Wannier function, V is potential
complex(8), allocatable :: s_wvlm(:,:,:,:)
!----------!
! energies !
!----------!
! {n,n',T} list of Wannier transitions used to compute matrix elements
type(wannier_transitions) :: sic_wantran
! matrix elements of Wannier potential <W_{n0}|V_{n0}|W_{n'T}>
complex(8), allocatable :: sic_vme(:)
! cutoff distance for SIC marix elements <W_n|V_n|W_{n'T}>
real(8) sic_me_cutoff
data sic_me_cutoff/0.1d0/
! LDA Hamiltonian in k-space in Wannier basis 
complex(8), allocatable :: sic_wan_h0k(:,:,:)
! LDA energies of Wannier functions
real(8), allocatable :: sic_wan_e0(:)
! DFT energy before the SIC correction
real(8) engytot0
! total energy correction
real(8) sic_energy_tot
data sic_energy_tot/0.d0/
! potential contribution
real(8) sic_energy_pot
data sic_energy_pot/0.d0/
! kinetic contribution
real(8) sic_energy_kin
data sic_energy_kin/0.d0/
!--------------!
! LAPW liaison !
!--------------!
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
! dot-product <W_{n\sigma}|f_{jk}>, f_{jk} is the first-variational Bloch state 
complex(8), allocatable :: sic_wb(:,:,:,:)
! dot-product <(W*V)_{n\sigma}|f_{jk}>, f_{jk} is the first-variational Bloch state 
complex(8), allocatable :: sic_wvb(:,:,:,:)
! Bloch-sum of Wannier functions in muffin-tins
complex(8), allocatable :: s_wkmt(:,:,:,:,:,:)
! Bloch-sum of Wannier functions in the interstitial expanded in plane waves
complex(8), allocatable :: s_wkit(:,:,:,:)
! Bloch-sum of W*V
complex(8), allocatable :: s_wvkmt(:,:,:,:,:,:)
! Bloch-sum of W*V in the interstitial expanded in plane waves
complex(8), allocatable :: s_wvkit(:,:,:,:)


! number of SIC iterations
integer nsclsic
data nsclsic/3/
! current SIC iteration
integer isclsic
data isclsic/0/
integer, allocatable :: sic_apply(:)
integer, allocatable :: sicw(:,:)
logical :: tsic_wv
data tsic_wv/.false./
integer sic_debug_level
data sic_debug_level/0/

! maximum number of translation vectors
integer, parameter :: sic_maxvtl=1000

type t_sic_blochsum
! total number of translations
  integer ntr
! translation vectors in lattice coordinates
  integer, allocatable :: vtl(:,:)
! translation vectors in Cartesian coordinates
  real(8), allocatable :: vtc(:,:)
end type t_sic_blochsum

type(t_sic_blochsum) :: sic_blochsum

integer, parameter :: nwanprop=14
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
integer, parameter :: wp_spread_x=11
integer, parameter :: wp_spread_y=12
integer, parameter :: wp_spread_z=13
integer, parameter :: wp_normrho=14

contains

subroutine s_get_wanval(x,wanval,itp,ir,rmax,rcutoff)
use modmain
use mod_nrkp
use mod_wannier
use mod_madness
implicit none
real(8), intent(inout) :: x(3)
complex(8), intent(out) :: wanval(nspinor)
integer, optional, intent(in) :: itp
integer, optional, intent(in) :: ir
real(8), optional, intent(in) :: rmax
real(8), optional, intent(in) :: rcutoff
!
integer is,ia,ias,ir0,io,l,j,i,lm,ig,ispn
integer ntr(3),ik,ikloc,ic
real(8) x0(3),vtc(3),vr0(3),r0,tp(2),t1
real(8) ur(0:lmaxvr,nufrmax),dr
complex(8) zt1,zt2,ylm(lmmaxvr)
logical, external :: vrinmt2
complex(8), external :: ylm_val
complex(8) expigr(m_ngvec)
complex(8) zm1(lmmaxvr,nufrmax,nspinor),zm2(nspinor)
!
wanval=zzero
if (present(itp).and.present(ir)) then
  x0(:)=s_x(:,itp)*s_r(ir)
  x(:)=x0(:)+m_wanpos(:)
else
  x0(:)=x(:)-m_wanpos(:)
endif
if (present(rmax)) then
  if (sum(x0(:)**2).gt.(rmax**2+1d-8)) then
    write(*,'("Error(s_get_wanval2): outside of rmax")')
    write(*,'("  rmax : ",G18.10)')rmax
    write(*,'("   x   : ",3G18.10)')x0
    write(*,'("  |x|  : ",G18.10)')sqrt(sum(x0(:)**2))
    call pstop
  endif
endif
if (present(rcutoff)) then
  if (sum(x0(:)**2).gt.(rcutoff**2)) return
endif

if (vrinmt2(x,is,ia,ntr,ir0,vr0,dr)) then
  ias=idxas(ia,is)
  ic=ias2ic(ias)
  call sphcrd(vr0,r0,tp)
  call genylm(lmaxvr,tp,ylm)
  vtc(:)=ntr(1)*avec(:,1)+ntr(2)*avec(:,2)+ntr(3)*avec(:,3)
  ur=0.d0
  do l=0,lmaxvr
    do io=1,nufr(l,is)
      ur(l,io)=ufr(ir0,l,io,ic)+dr*(ufr(ir0+1,l,io,ic)-ufr(ir0,l,io,ic))
    enddo !io
  enddo !l
  zm1=zzero
  do ik=1,nkptnr
    zt1=exp(zi*dot_product(vkcnr(:,ik),vtc(:)))*wkptnr(ik)
    do ispn=1,nspinor
      zm1(:,:,ispn)=zm1(:,:,ispn)+zt1*m_wann_unkmt(:,:,ias,ispn,ik)
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
  do ig=1,m_ngvec
    expigr(ig)=exp(zi*dot_product(x(:),vgc(:,ig)))  
  enddo
  do ik=1,nkptnr
    zt1=wkptnr(ik)*exp(zi*dot_product(vkcnr(:,ik),x(:)))/sqrt(omega)
    zm2=zzero
! TODO: check if this can be optimized by switching order of indexes
    do ispn=1,nspinor
      do ig=1,m_ngknr(ik)
      zt2=expigr(m_igkignr(ig,ik))
      !do ispn=1,nspinor
        zm2(ispn)=zm2(ispn)+zt2*m_wann_unkit(ig,ispn,ik)
      enddo
    enddo
    do ispn=1,nspinor
      wanval(ispn)=wanval(ispn)+zt1*zm2(ispn)
    enddo
  enddo
endif
return
end subroutine

subroutine s_spinor_func_val(x,f1lm,zval1,f2lm,zval2,rmax,rcutoff)
use modmain
implicit none
! arguments
real(8), intent(in) :: x(3)
complex(8), intent(in) :: f1lm(lmmaxwan,s_nr,nspinor)
complex(8), intent(out) :: zval1(nspinor)
complex(8), optional, intent(in) :: f2lm(lmmaxwan,s_nr,nspinor)
complex(8), optional, intent(out) :: zval2(nspinor)
real(8), optional, intent(in) :: rmax
real(8), optional, intent(in) :: rcutoff
! local variables
integer ir1,lm,ir,ispn
real (8) x0,tp(2),dx
complex(8) ylm(lmmaxwan)
complex(8) z1,z2
!
if (present(rmax)) then
  if (sum(x(:)**2).gt.(rmax**2+1d-8)) then
    write(*,'("Error(s_spinor_func_val): outside of rmax")')
    write(*,'("  rmax : ",G18.10)')rmax
    write(*,'("   x   : ",3G18.10)')x
    write(*,'("  |x|  : ",G18.10)')sqrt(sum(x(:)**2))
    call pstop
  endif
endif
if (present(rcutoff)) then
  if (sum(x(:)**2).gt.(rcutoff**2)) then
    zval1=zzero
    if (present(zval2)) zval2=zzero
    return
  endif
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
  do lm=1,lmmaxwan
    z1=z1+(f1lm(lm,ir1,ispn)+dx*(f1lm(lm,ir1+1,ispn)-f1lm(lm,ir1,ispn)))*ylm(lm)
  enddo
  zval1(ispn)=z1
  if (present(zval2)) then
    z2=zzero
    do lm=1,lmmaxwan
      z2=z2+(f2lm(lm,ir1,ispn)+dx*(f2lm(lm,ir1+1,ispn)-f2lm(lm,ir1,ispn)))*ylm(lm)
    enddo
    zval2(ispn)=z2
  endif
enddo
return
end subroutine 

complex(8) function s_spinor_dotp(pos1,pos2,f1lm,f2lm)
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
    do ir=1,s_nr_min
      zprod=zprod+zdotc(lmmaxwan,f1lm(1,ir,ispn),1,f2lm(1,ir,ispn),1)*s_rw(ir)
    enddo
  enddo
else
  allocate(f1tp_(s_ntp,s_nr))
  allocate(f2tp_(s_ntp,s_nr,nspinor))
  !allocate(f2lm_(lmmaxwan,s_nr))
  f2tp_=zzero
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(itp,x1,x2)
  do ir=1,s_nr_min
    do itp=1,s_ntp
      x1(:)=s_x(:,itp)*s_r(ir)
      x2(:)=pos1(:)+x1(:)-pos2(:)
      call s_spinor_func_val(x2,f2lm,f2tp_(itp,ir,:),rcutoff=s_rmin)
    enddo
  enddo !ir
!$OMP END PARALLEL DO
  do ispn=1,nspinor
! convert to spherical coordinates
    call zgemm('T','N',s_ntp,s_nr,lmmaxwan,zone,s_ylmf,lmmaxwan,&
      f1lm(1,1,ispn),lmmaxwan,zzero,f1tp_,s_ntp)
    do ir=1,s_nr_min
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
s_spinor_dotp=zprod
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
real(8) t1,x(3),x2
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
allocate(rhotp(s_ntp,s_nr_min,nspinor))
allocate(rholm(lmmaxwan,s_nr_min,nspinor))
allocate(totrholm(lmmaxwan,s_nr_min))
allocate(vhalm(lmmaxwan,s_nr_min))
allocate(extp(s_ntp,s_nr_min))
allocate(ectp(s_ntp,s_nr_min))
allocate(exclm(lmmaxwan,s_nr_min))
allocate(vxtp(s_ntp,s_nr_min,nspinor))
allocate(vctp(s_ntp,s_nr_min,nspinor))
allocate(vxclm(lmmaxwan,s_nr_min,nspinor))

lmmaxwanloc=mpi_grid_map2(lmmaxwan,dims=(/dim_k,dim2/))

! TODO: memory optimizaton?
totrholm=0.d0
rholm=0.d0
do ispn=1,nspinor
  rhotp(:,1:s_nr_min,ispn)=abs(wantp(:,1:s_nr_min,ispn))**2
enddo
! charge density
! w(r)=\sum_{L} w_{L}(r) Y_{L}(t,p)
! rho(r) = \sum_{L2} rho_{L2}(r) R_{L2}(t,p) = 
!  = \sum_{L1,L3} w_{L1}^{*}(r) Y_{L1}^{*}(t,p) * w_{L3}(r) Y_{L3} (t,p)
! rho_{L2}(r) = \sum_{L1,L3} w_{L1}^{*}(r) w_{L3}(r)   <Y_{L1} |R_{L2}| Y_{L3}>
allocate(f1(s_nr_min,lmmaxwan,nspinor))
allocate(f2(s_nr_min,lmmaxwan,nspinor))
! rearrange wanlm in memory
do ispn=1,nspinor
  do lm=1,lmmaxwan
    f1(:,lm,ispn)=wanlm(lm,1:s_nr_min,ispn)
  enddo
enddo
f2=0.d0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(lm1,lm2,lm3,zt1,ispn)
do lmloc=1,lmmaxwanloc
  lm2=mpi_grid_map2(lmmaxwan,dims=(/dim_k,dim2/),loc=lmloc)
  do lm1=1,lmmaxwan
    do lm3=1,lmmaxwan
      zt1=gauntyry(lm2l(lm1),lm2l(lm2),lm2l(lm3),&
                   lm2m(lm1),lm2m(lm2),lm2m(lm3))
      if (abs(zt1).gt.1d-12) then
        do ispn=1,nspinor
          f2(:,lm2,ispn)=f2(:,lm2,ispn)+&
            dreal(dconjg(f1(:,lm1,ispn))*f1(:,lm3,ispn)*zt1)
        enddo !ispn
      endif
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
call mpi_grid_reduce(f2(1,1,1),s_nr_min*lmmaxwan*nspinor,all=.true.)
do ispn=1,nspinor
  do lm=1,lmmaxwan
    rholm(lm,:,ispn)=f2(:,lm,ispn)
  enddo
  totrholm(:,:)=totrholm(:,:)+rholm(:,:,ispn)
enddo
deallocate(f1,f2)
! norm of total charge density
t1=0.d0
do ir=1,s_nr_min
  t1=t1+totrholm(1,ir)*s_rw(ir)
enddo
wanprop(wp_normrho)=t1*fourpi*y00
! estimate the quadratic spread <r^2>-<r>^2
x2=0.d0
x=0.d0
! Ry=-\frac{1}{2} \sqrt{\frac{3}{\pi }} \sin (t) \sin (p)
! Rz= \frac{1}{2} \sqrt{\frac{3}{\pi }} \cos (t)
! Rx=-\frac{1}{2} \sqrt{\frac{3}{\pi }} \sin (t) \cos (p) 
do ir=1,s_nr_min
  x(1)=x(1)-2.d0*s_r(ir)*sqrt(pi/3.d0)*totrholm(4,ir)*s_rw(ir)
  x(2)=x(2)-2.d0*s_r(ir)*sqrt(pi/3.d0)*totrholm(2,ir)*s_rw(ir)
  x(3)=x(3)+2.d0*s_r(ir)*sqrt(pi/3.d0)*totrholm(3,ir)*s_rw(ir)
  x2=x2+2.d0*(s_r(ir)**2)*sqrt(pi)*totrholm(1,ir)*s_rw(ir)
enddo
wanprop(wp_spread)=x2-dot_product(x,x)
wanprop(wp_spread_x)=x(1)
wanprop(wp_spread_y)=x(2)
wanprop(wp_spread_z)=x(3)
! compute Hartree potential
vhalm=0.d0
do lmloc=1,lmmaxwanloc
  lm=mpi_grid_map2(lmmaxwan,dims=(/dim_k,dim2/),loc=lmloc)
  l=lm2l(lm)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jr,t1)
  do ir=1,s_nr_min
    t1=0.d0
    do jr=1,ir
      t1=t1+(s_r(jr)**l/s_r(ir)**(l+1))*totrholm(lm,jr)*s_rw(jr)
    enddo
    do jr=ir+1,s_nr_min
      t1=t1+(s_r(ir)**l/s_r(jr)**(l+1))*totrholm(lm,jr)*s_rw(jr)
    enddo
    vhalm(lm,ir)=t1*fourpi/(2*l+1)
  enddo
!$OMP END PARALLEL DO
enddo
call mpi_grid_reduce(vhalm(1,1),lmmaxwan*s_nr_min,all=.true.)
! compute XC potential
if (spinpol) then
  call xcifc(xctype,n=s_ntp*s_nr_min,rhoup=rhotp(1,1,1),rhodn=rhotp(1,1,2),&
    ex=extp,ec=ectp,vxup=vxtp(1,1,1),vxdn=vxtp(1,1,2),vcup=vctp(1,1,1),&
    vcdn=vctp(1,1,2))
else
  call xcifc(xctype,n=s_ntp*s_nr_min,rho=rhotp,ex=extp,ec=ectp,vx=vxtp,vc=vctp)
endif
! save XC energy density in extp
extp(:,:)=extp(:,:)+ectp(:,:)
! expand in real spherical harmonics
call sic_rbsht(s_nr_min,extp,exclm) 
! save XC potential in vxtp and expand in real spherical harmonics   
do ispn=1,nspinor
  vxtp(:,:,ispn)=vxtp(:,:,ispn)+vctp(:,:,ispn)
  call sic_rbsht(s_nr_min,vxtp(1,1,ispn),vxclm(1,1,ispn))
enddo
! compute vha=<V_h|rho>
wanprop(wp_vha)=0.d0
do ir=1,s_nr_min
  wanprop(wp_vha)=wanprop(wp_vha)+&
    ddot(lmmaxwan,totrholm(1,ir),1,vhalm(1,ir),1)*s_rw(ir)
enddo
! compute exc=<E_xc|rho>
wanprop(wp_exc)=0.d0
do ir=1,s_nr_min
  wanprop(wp_exc)=wanprop(wp_exc)+&
    ddot(lmmaxwan,totrholm(1,ir),1,exclm(1,ir),1)*s_rw(ir)
enddo
! compute vxc=<W_n|V_xc|W_n>; in the collinear case this is 
!  \sum_{\sigma} <V_xc^{\sigma}|rho${\sigma}>
wanprop(wp_vxc)=0.d0
do ispn=1,nspinor
  do ir=1,s_nr_min
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
allocate(f1(s_nr_min,lmmaxwan,nspinor))
allocate(f2(s_nr_min,lmmaxwan,nspinor))
allocate(f3(s_nr_min,lmmaxwan,nspinor))
! rearrange arrays in memory
do ispn=1,nspinor
  do lm=1,lmmaxwan
    f1(:,lm,ispn)=wanlm(lm,1:s_nr_min,ispn)
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
          f3(:,lm3,ispn)=f3(:,lm3,ispn)-f1(:,lm1,ispn)*f2(:,lm2,ispn)*zt1
        enddo !ispn
      endif
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
call mpi_grid_reduce(f3(1,1,1),s_nr_min*lmmaxwan*nspinor,all=.true.) 
wvlm=zzero
do ispn=1,nspinor
  do lm=1,lmmaxwan
    wvlm(lm,1:s_nr_min,ispn)=f3(1:s_nr_min,lm,ispn)
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
  zsumir=zdotc(ng,zfir1,1,zfir2,1)
  !do ig=1,ng
  !  zsumir=zsumir+dconjg(zfir1(ig))*zfir2(ig)
  !enddo
else
  do ir=1,ngrtot
    zsumir=zsumir+cfunir(ir)*dconjg(zfir1(ir))*zfir2(ir)
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

subroutine sic_zbsht(nr,zftp,zflm)
use modmain
implicit none
integer, intent(in) :: nr
complex(8), intent(in) :: zftp(s_ntp,nr)
complex(8), intent(inout) :: zflm(lmmaxwan,nr)
!
integer iter,nrloc,roffs
real(8), allocatable :: tdiff(:)
complex(8), allocatable :: zftp1(:,:)
!
nrloc=mpi_grid_map(nr,dim_k,offs=roffs)
zflm=zzero
! convert to spherical harmonics
call zgemm('T','N',lmmaxwan,nrloc,s_ntp,zone,s_ylmb,s_ntp,&
  zftp(1,roffs+1),s_ntp,zzero,zflm(1,roffs+1),lmmaxwan)
allocate(zftp1(s_ntp,nrloc))
allocate(tdiff(sic_bsht_niter))
tdiff=0.d0
do iter=1,sic_bsht_niter
! convert back to spherical coordinates
  call zgemm('T','N',s_ntp,nrloc,lmmaxwan,zone,s_ylmf,lmmaxwan,zflm(1,1+roffs),&
    lmmaxwan,zzero,zftp1,s_ntp)
  zftp1(:,1:nrloc)=zftp(:,roffs+1:roffs+nrloc)-zftp1(:,1:nrloc)
  tdiff(iter)=sum(abs(zftp1))
! add difference to spherical harmonic expanson
  call zgemm('T','N',lmmaxwan,nrloc,s_ntp,zone,s_ylmb,s_ntp,&
    zftp1,s_ntp,zone,zflm(1,roffs+1),lmmaxwan)
enddo
deallocate(zftp1)
if (sic_bsht_niter.gt.0) then
  call mpi_grid_reduce(tdiff(1),sic_bsht_niter,dims=(/dim_k/))
endif
call mpi_grid_reduce(zflm(1,1),lmmaxwan*nr,dims=(/dim_k/),all=.true.)
if (mpi_grid_root().and.sic_bsht_niter.gt.1) then
  if (tdiff(sic_bsht_niter).gt.tdiff(sic_bsht_niter-1).and.&
      tdiff(sic_bsht_niter).gt.1d-6) then
    write(*,'("Warning(sic_zbsht): difference of functions at each iteration")')
    write(*,'(255F12.6)')tdiff
  endif
endif
deallocate(tdiff)
return
end subroutine

subroutine sic_rbsht(nr,ftp,flm)
use modmain
implicit none
integer, intent(in) :: nr
real(8), intent(in) :: ftp(s_ntp,nr)
real(8), intent(inout) :: flm(lmmaxwan,nr)
!
integer iter,nrloc,roffs
real(8), allocatable :: tdiff(:)
real(8), allocatable :: ftp1(:,:)
!
nrloc=mpi_grid_map(nr,dim_k,offs=roffs)
flm=0.d0
! convert to spherical harmonics
call dgemm('T','N',lmmaxwan,nrloc,s_ntp,1.d0,s_rlmb,s_ntp,&
  ftp(1,roffs+1),s_ntp,0.d0,flm(1,roffs+1),lmmaxwan)
allocate(ftp1(s_ntp,nrloc))
allocate(tdiff(sic_bsht_niter))
tdiff=0.d0
do iter=1,sic_bsht_niter
! convert back to spherical coordinates
  call dgemm('T','N',s_ntp,nrloc,lmmaxwan,1.d0,s_rlmf,lmmaxwan,flm(1,1+roffs),&
    lmmaxwan,0.d0,ftp1,s_ntp)
  ftp1(:,1:nrloc)=ftp(:,roffs+1:roffs+nrloc)-ftp1(:,1:nrloc)
  tdiff(iter)=sum(abs(ftp1))
! add spherical harmonic expanson of the difference to the total expansion
  call dgemm('T','N',lmmaxwan,nrloc,s_ntp,1.d0,s_rlmb,s_ntp,&
    ftp1,s_ntp,1.d0,flm(1,roffs+1),lmmaxwan)
enddo
deallocate(ftp1)
if (sic_bsht_niter.gt.0) then
  call mpi_grid_reduce(tdiff(1),sic_bsht_niter,dims=(/dim_k/))
endif
call mpi_grid_reduce(flm(1,1),lmmaxwan*nr,dims=(/dim_k/),all=.true.)
if (mpi_grid_root().and.sic_bsht_niter.gt.1) then
  if (tdiff(sic_bsht_niter).gt.tdiff(sic_bsht_niter-1).and.&
      tdiff(sic_bsht_niter).gt.1d-6) then
    write(*,'("Warning(sic_rbsht): difference of functions at each iteration")')
    write(*,'(255F12.6)')tdiff
  endif
endif
deallocate(tdiff)
return
end subroutine

end module
