module mod_sic
use mod_wannier
!---------------------!
! spherical expansion !
!---------------------!
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
! number of radial points for s_rmin
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
! dot-product <W_{n\sigma}|b_{jk}>, b_{jk} is the basis Bloch orbital 
complex(8), allocatable :: sic_wb(:,:,:,:)
! dot-product <(W*V)_{n\sigma}|b_{jk}>, b_{jk} is the basis Bloch orbital
complex(8), allocatable :: sic_wvb(:,:,:,:)
! lm-expanded Bloch-sum of Wannier functions in muffin-tins
! note: order of indexes is optimized for dot-product with (L)APW basis
complex(8), allocatable :: s_wkmt(:,:,:,:,:,:)
! Bloch-sum of Wannier functions in the interstitial expanded in plane waves
complex(8), allocatable :: s_wkit(:,:,:,:)
! lm-expanded Bloch-sum of W*V
complex(8), allocatable :: s_wvkmt(:,:,:,:,:,:)
! Bloch-sum of W*V in the interstitial expanded in plane waves
complex(8), allocatable :: s_wvkit(:,:,:,:)
! number of translations for Bloch-sums
integer sic_ntr
! translation vectors in lattice coordinates
integer, allocatable :: sic_vtl(:,:)
! translation vectors in Cartesian coordinates
real(8), allocatable :: sic_vtc(:,:)
! unitary matrix for the localization criterion
complex(8), allocatable :: sic_wan_umtrx(:,:,:)
real(8) sic_umtrx_eps
data sic_umtrx_eps/1.d0/
real(8) sic_u0_eps
data sic_u0_eps/1.d0/
integer sic_niter_umtrx
data sic_niter_umtrx/1/
integer sic_niter_u0
data sic_niter_u0/1/

logical tsicsv
data tsicsv/.false./

! .true. if sic branch is activated 
logical sic
data sic/.false./
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
integer is,ia,ias,ir0,io,l,lm,ig,ispn
integer ntr(3),ik,ic
real(8) x0(3),vtc(3),vr0(3),r0,tp(2)
real(8) ur(0:lmaxapw,nufrmax),dr
complex(8) zt1,zt2,ylm(lmmaxapw)
logical, external :: vrinmt2
complex(8), external :: ylm_val
complex(8) expigr(m_ngvec)
complex(8) zm1(lmmaxapw,nufrmax,nspinor),zm2(nspinor)
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
  call genylm(lmaxapw,tp,ylm)
  vtc(:)=ntr(1)*avec(:,1)+ntr(2)*avec(:,2)+ntr(3)*avec(:,3)
  ur=0.d0
  do l=0,lmaxapw
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
    do lm=1,lmmaxapw
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
complex(8), optional, intent(out) :: zfrac(2)
!
integer is,ias,ir,itp
complex(8) zsumir,zsummt,zt1
complex(8), external :: zdotc
! interstitial contribution
zsumir=zzero
if (tpw) then
  zsumir=zdotc(ng,zfir1,1,zfir2,1)
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
  zfrac(1)=zsummt
  zfrac(2)=zsumir
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
allocate(tdiff(sic_bsht_niter))
tdiff=0.d0
if (nrloc.gt.0) then
! convert to spherical harmonics
  call zgemm('T','N',lmmaxwan,nrloc,s_ntp,zone,s_ylmb,s_ntp,&
    zftp(1,roffs+1),s_ntp,zzero,zflm(1,roffs+1),lmmaxwan)
  allocate(zftp1(s_ntp,nrloc))
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
endif
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
allocate(tdiff(sic_bsht_niter))
tdiff=0.d0
if (nrloc.gt.0) then
! convert to spherical harmonics
  call dgemm('T','N',lmmaxwan,nrloc,s_ntp,1.d0,s_rlmb,s_ntp,&
    ftp(1,roffs+1),s_ntp,0.d0,flm(1,roffs+1),lmmaxwan)
  allocate(ftp1(s_ntp,nrloc))
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
endif
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


subroutine sic_genbprj(ikloc,evecfv,apwalm,wfsvmt,wfsvit)
use modmain
implicit none
! arguments
integer, intent(in) :: ikloc
complex(8), optional, intent(in) :: evecfv(nmatmax,nstfv)
complex(8), optional, intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), optional, intent(in) :: wfsvmt(lmmaxapw,nufrmax,natmtot,nspinor,nstsv)
complex(8), optional, intent(in) :: wfsvit(ngkmax,nspinor,nstsv)
! local variables
integer ik,i,j,ispn,l,m,lm,io,ilo,ias,is,ic,ist,ir,ig
complex(8) z1,z2
complex(8), allocatable :: wfmt1(:,:,:)
complex(8), allocatable :: wfmt2(:,:,:)
complex(8), allocatable :: zv1(:),zv2(:)
!
sic_wb(:,:,:,ikloc)=zzero
sic_wvb(:,:,:,ikloc)=zzero
if (.not.tsic_wv) return 
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
! project to first-variational states
if (tsveqn)  then
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(wfmt1,wfmt2,ias,j,ispn)
  allocate(wfmt1(lmmaxapw,nrmtmax,natmtot))
  allocate(wfmt2(lmmaxapw,nrmtmax,natmtot))
!$OMP DO
  do ist=1,nstfv
    wfmt1=zzero
! generate first-variational wave function
    do ias=1,natmtot
      call wavefmt(1,lmaxapw,ias2is(ias),ias2ia(ias),ngk(1,ik),apwalm,&
        evecfv(1,ist),lmmaxapw,wfmt1(1,1,ias))
    enddo
    do j=1,sic_wantran%nwan
      do ispn=1,nspinor
! rearrange indexes
        do ias=1,natmtot
          do lm=1,lmmaxapw
            wfmt2(lm,:,ias)=s_wkmt(:,lm,ias,ispn,j,ikloc)
          enddo
        enddo
        sic_wb(j,ist,ispn,ikloc)=s_zfinp(.true.,.true.,lmmaxapw,ngk(1,ik),&
          wfmt2,wfmt1,s_wkit(1,ispn,j,ikloc),evecfv(1,ist))
! rearrange indexes
        do ias=1,natmtot
          do lm=1,lmmaxapw
            wfmt2(lm,:,ias)=s_wvkmt(:,lm,ias,ispn,j,ikloc)
          enddo
        enddo
        sic_wvb(j,ist,ispn,ikloc)=s_zfinp(.true.,.true.,lmmaxapw,ngk(1,ik),&
          wfmt2,wfmt1,s_wvkit(1,ispn,j,ikloc),evecfv(1,ist))
      enddo
    enddo
  enddo !ist
!$OMP END DO
  deallocate(wfmt1,wfmt2)
!$OMP END PARALLEL
else
  if (tsicsv) then
! project to second-variational states
    do j=1,sic_wantran%nwan
      do ispn=1,nspinor
        do ias=1,natmtot
          is=ias2is(ias)
          ic=ias2ic(ias)
          do lm=1,lmmaxapw
            l=lm2l(lm)
            do io=1,nufr(l,is)
              z1=zzero
              z2=zzero
              do ir=1,nrmt(is)
                z1=z1+dconjg(s_wkmt(ir,lm,ias,ispn,j,ikloc))*ufr(ir,l,io,ic)*mt_rw(ir,is)
                z2=z2+dconjg(s_wvkmt(ir,lm,ias,ispn,j,ikloc))*ufr(ir,l,io,ic)*mt_rw(ir,is)
              enddo !ir
              do ist=1,nstsv
                sic_wb(j,ist,1,ikloc)=sic_wb(j,ist,1,ikloc)+z1*wfsvmt(lm,io,ias,ispn,ist)
                sic_wvb(j,ist,1,ikloc)=sic_wvb(j,ist,1,ikloc)+z2*wfsvmt(lm,io,ias,ispn,ist)
              enddo !ist
            enddo !io
          enddo !lm
        enddo !ias
        do ist=1,nstsv
          do ig=1,ngk(1,ik)
            sic_wb(j,ist,1,ikloc)=sic_wb(j,ist,1,ikloc)+&
              dconjg(s_wkit(ig,ispn,j,ikloc))*wfsvit(ig,ispn,ist)
            sic_wvb(j,ist,1,ikloc)=sic_wvb(j,ist,1,ikloc)+&
              dconjg(s_wvkit(ig,ispn,j,ikloc))*wfsvit(ig,ispn,ist)
          enddo !ig
        enddo !ist
      enddo !ispn
    enddo !j  
  else
! project to (L)APW basis
    allocate(zv1(ngk(1,ik)))
    allocate(zv2(ngk(1,ik)))
    do j=1,sic_wantran%nwan
      do ispn=1,nspinor
! interstitial contribution from APW
        do ig=1,ngk(1,ik)
          zv1(ig)=dconjg(s_wkit(ig,ispn,j,ikloc))
          zv2(ig)=dconjg(s_wvkit(ig,ispn,j,ikloc))
        enddo
! muffin-tin contribution from APW
        do ias=1,natmtot
          is=ias2is(ias)
          do lm=1,lmmaxapw
            l=lm2l(lm)
            do io=1,apword(l,is)
              z1=zzero
              z2=zzero
              do ir=1,nrmt(is)
                z1=z1+dconjg(s_wkmt(ir,lm,ias,ispn,j,ikloc))*&
                  apwfr(ir,1,io,l,ias)*mt_rw(ir,is)
                z2=z2+dconjg(s_wvkmt(ir,lm,ias,ispn,j,ikloc))*&
                  apwfr(ir,1,io,l,ias)*mt_rw(ir,is)
              enddo !ir
              do ig=1,ngk(1,ik)
                zv1(ig)=zv1(ig)+z1*apwalm(ig,io,lm,ias)
                zv2(ig)=zv2(ig)+z2*apwalm(ig,io,lm,ias)
              enddo !ig
            enddo !io
          enddo !lm
! muffin-tin contribution from l.o.
          do ilo=1,nlorb(is)
            l=lorbl(ilo,is)
            do m=-l,l
              lm=idxlm(l,m)
              i=ngk(1,ik)+idxlo(lm,ilo,ias)
              z1=zzero
              z2=zzero
              do ir=1,nrmt(is)
                z1=z1+dconjg(s_wkmt(ir,lm,ias,ispn,j,ikloc))*&
                  lofr(ir,1,ilo,ias)*mt_rw(ir,is)
                z2=z2+dconjg(s_wvkmt(ir,lm,ias,ispn,j,ikloc))*&
                  lofr(ir,1,ilo,ias)*mt_rw(ir,is)
              enddo
              sic_wb(j,i,ispn,ikloc)=z1
              sic_wvb(j,i,ispn,ikloc)=z2
            enddo !m
          enddo !ilo
        enddo !ias
        sic_wb(j,1:ngk(1,ik),ispn,ikloc)=zv1(1:ngk(1,ik))
        sic_wvb(j,1:ngk(1,ik),ispn,ikloc)=zv2(1:ngk(1,ik))
      enddo !ispn
    enddo !j
    deallocate(zv1,zv2)
  endif
endif
return
end subroutine

end module
