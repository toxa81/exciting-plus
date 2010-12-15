module mod_sic
use mod_wannier

! matrix elements of Wannier potential <W_{n0}|V_{n0}|W_{n'T}>
complex(8), allocatable :: vwanme(:)
! LDA Hamiltonian in k-space in Wannier basis 
complex(8), allocatable :: sic_wann_h0k(:,:,:)
! LDA energies of Wannier functions
real(8), allocatable :: sic_wann_e0(:)
integer :: nsclsic
data nsclsic/3/
integer :: isclsic
data isclsic/0/
real(8) :: etot0 
real(8) :: sic_etot_correction
data sic_etot_correction/0.d0/
real(8) :: sic_wan_cutoff
data sic_wan_cutoff/6.d0/
real(8) :: sic_me_cutoff
data sic_me_cutoff/0.1d0/
complex(8), allocatable :: sic_wb(:,:,:,:)
complex(8), allocatable :: sic_wvb(:,:,:,:)

integer, allocatable :: sic_apply(:)
integer, allocatable :: sicw(:,:)

logical :: tsic_wv
data tsic_wv/.false./
logical :: tsic_arrays_allocated
data tsic_arrays_allocated/.false./

! maximum number of translation vectors
integer, parameter :: sic_maxvtl=1000

integer :: ngrloc
integer :: ngrlocmax
integer :: groffs
integer :: nmtloc
integer :: nmtlocmax
integer :: mtoffs
! weights for radial integration
real(8), allocatable :: rmtwt(:)

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
! Wannier functions
  complex(8), allocatable :: wanmt(:,:,:,:,:)
  complex(8), allocatable :: wanir(:,:,:,:)
! product of a Wannier function with it's potential
  complex(8), allocatable :: wvmt(:,:,:,:,:)
  complex(8), allocatable :: wvir(:,:,:,:)
! .true. if Wannier function is expanded inside muffin-tin in the given cell
  logical, allocatable :: twanmt(:,:,:)
! .true. if at least one Wannier function is expanded inside muffin-tin in 
!   the given unit cell
  logical, allocatable :: twanmtuc(:,:)
end type t_sic_orbitals

type(t_sic_orbitals) :: sic_orbitals


!interface sic_copy_mt
!  module procedure sic_copy_mt_z,sic_copy_mt_d
!end interface
!
!interface sic_copy_ir
!  module procedure sic_copy_ir_z
!end interface

contains

subroutine sic_copy_mt_z(tfrwrd,ld,fmt1,fmt2)
use modmain
implicit none
logical, intent(in) :: tfrwrd
integer, intent(in) :: ld
complex(8), intent(inout) :: fmt1(ld,nrmtmax*natmtot)
complex(8), intent(inout) :: fmt2(ld,nmtloc)
if (tfrwrd) then
  fmt2(1:ld,1:nmtloc)=fmt1(1:ld,mtoffs+1:mtoffs+nmtloc)
else
  fmt1(1:ld,mtoffs+1:mtoffs+nmtloc)=fmt2(1:ld,1:nmtloc)
endif
return
end subroutine

subroutine sic_copy_mt_d(tfrwrd,ld,fmt1,fmt2)
use modmain
implicit none
logical, intent(in) :: tfrwrd
integer, intent(in) :: ld
real(8), intent(inout) :: fmt1(ld,nrmtmax*natmtot)
real(8), intent(inout) :: fmt2(ld,nmtloc)
if (tfrwrd) then
  fmt2(1:ld,1:nmtloc)=fmt1(1:ld,mtoffs+1:mtoffs+nmtloc)
else
  fmt1(1:ld,mtoffs+1:mtoffs+nmtloc)=fmt2(1:ld,1:nmtloc)
endif
return
end subroutine

subroutine sic_copy_ir_z(tfrwrd,fir1,fir2)
use modmain
implicit none
logical, intent(in) :: tfrwrd
complex(8), intent(inout) :: fir1(ngrtot)
complex(8), intent(inout) :: fir2(ngrloc)
if (tfrwrd) then
  fir2(1:ngrloc)=fir1(groffs+1:groffs+ngrloc)
else
  fir1(groffs+1:groffs+ngrloc)=fir2(1:ngrloc)
endif
return
end subroutine

subroutine sic_copy_ir_d(tfrwrd,fir1,fir2)
use modmain
implicit none
logical, intent(in) :: tfrwrd
real(8), intent(inout) :: fir1(ngrtot)
real(8), intent(inout) :: fir2(ngrloc)
if (tfrwrd) then
  fir2(1:ngrloc)=fir1(groffs+1:groffs+ngrloc)
else
  fir1(groffs+1:groffs+ngrloc)=fir2(1:ngrloc)
endif
return
end subroutine



subroutine sic_copy_mt_z_2(ld,nmtloc1,mtoffs1,fmt1,fmt2)
use modmain
implicit none
integer, intent(in) :: ld
integer, intent(in) :: nmtloc1
integer, intent(in) :: mtoffs1
complex(8), intent(in) :: fmt1(ld,nrmtmax*natmtot)
complex(8), intent(out) :: fmt2(ld,nmtloc1)
fmt2(1:ld,1:nmtloc1)=fmt1(1:ld,mtoffs1+1:mtoffs1+nmtloc1)
return
end subroutine

subroutine sic_copy_ir_z_2(ngrloc1,groffs1,fir1,fir2)
use modmain
implicit none
integer, intent(in) :: ngrloc1
integer, intent(in) :: groffs1
complex(8), intent(in) :: fir1(ngrtot)
complex(8), intent(inout) :: fir2(ngrloc1)
fir2(1:ngrloc1)=fir1(groffs1+1:groffs1+ngrloc1)
return
end subroutine



! compute <f1_0|f2_T>=\int_{-\inf}^{\inf} f1^{*}(r)f2(r-T)dr = 
!   = \sum_{R} \int_{\Omega} f1^{*}(r+R)f2(r+R-T)dr
complex(8) function sic_dot_ll(fmt1,fir1,fmt2,fir2,t,tfmt1,tfmt2)
use modmain
implicit none
! arguments
complex(8), intent(in) :: fmt1(lmmaxvr,nmtloc,sic_orbitals%ntr)
complex(8), intent(in) :: fir1(ngrloc,sic_orbitals%ntr)
complex(8), intent(in) :: fmt2(lmmaxvr,nmtloc,sic_orbitals%ntr)
complex(8), intent(in) :: fir2(ngrloc,sic_orbitals%ntr)
integer, intent(in) :: t(3)
logical, intent(in) :: tfmt1(natmtot,sic_orbitals%ntr)
logical, intent(in) :: tfmt2(natmtot,sic_orbitals%ntr)
! local variables
integer v1(3),v2(3),jt,ir,ias,it,i
complex(8) zdotmt,zdotir,zt1
complex(8), external :: zdotc
zdotmt=zzero
zdotir=zzero
do it=1,sic_orbitals%ntr
  v1(:)=sic_orbitals%vtl(:,it)
  v2(:)=v1(:)-t(:)
  if (v2(1).ge.sic_orbitals%tlim(1,1).and.v2(1).le.sic_orbitals%tlim(2,1).and.&
      v2(2).ge.sic_orbitals%tlim(1,2).and.v2(2).le.sic_orbitals%tlim(2,2).and.&
      v2(3).ge.sic_orbitals%tlim(1,3).and.v2(3).le.sic_orbitals%tlim(2,3)) then
    jt=sic_orbitals%ivtit(v2(1),v2(2),v2(3))
    if (jt.ne.-1) then
      do i=1,nmtloc
        ias=(mtoffs+i-1)/nrmtmax+1
        if (tfmt1(ias,it).and.tfmt2(ias,jt)) then
          zdotmt=zdotmt+rmtwt(i)*zdotc(lmmaxvr,fmt1(:,i,it),1,&
            fmt2(:,i,jt),1)
        endif
      enddo
      do ir=1,ngrloc
        zdotir=zdotir+cfunir(ir+groffs)*dconjg(fir1(ir,it))*fir2(ir,jt)
      enddo
    endif
  endif
enddo
zt1=zdotmt+zdotir*omega/dble(ngrtot)
call mpi_grid_reduce(zt1,all=.true.)
sic_dot_ll=zt1
return
end function

! compute <f1|f2|f3>=\int_{-\inf}^{\inf} f1^{*}(r)f2(r)f3(r)dr = 
!   = \sum_{R} \int_{\Omega} f1^{*}(r+R)f2(r+R)f3(r+R)dr
!  f1,f3 are complex
!  f2 is real
complex(8) function sic_int_zdz(f1mt,f1ir,f2mt,f2ir,f3mt,f3ir,tfmtuc)
use modmain
implicit none
! arguments
complex(8), intent(in) :: f1mt(lmmaxvr,nmtloc,sic_orbitals%ntr)
complex(8), intent(in) :: f1ir(ngrloc,sic_orbitals%ntr)
real(8), intent(in) :: f2mt(lmmaxvr,nmtloc,sic_orbitals%ntr)
real(8), intent(in) :: f2ir(ngrloc,sic_orbitals%ntr)
complex(8), intent(in) :: f3mt(lmmaxvr,nmtloc,sic_orbitals%ntr)
complex(8), intent(in) :: f3ir(ngrloc,sic_orbitals%ntr)
logical, intent(in) :: tfmtuc(sic_orbitals%ntr)
! local variables
complex(8), allocatable :: f1mt_(:,:),f3mt_(:,:)
real(8), allocatable :: f2mt_(:,:)
complex(8) zsummt,zsumir
integer it,ir,lm,lm1,lm2,lm3
complex(8) zt1
complex(8), external :: gauntyry

allocate(f1mt_(nmtloc,lmmaxvr))
allocate(f2mt_(nmtloc,lmmaxvr))
allocate(f3mt_(nmtloc,lmmaxvr))
zsummt=zzero
zsumir=zzero
do it=1,sic_orbitals%ntr
  if (tfmtuc(it)) then
! muffin-tin part
    do lm=1,lmmaxvr
      f1mt_(:,lm)=f1mt(lm,:,it)
      f2mt_(:,lm)=f2mt(lm,:,it)
      f3mt_(:,lm)=f3mt(lm,:,it)
    enddo
    do lm1=1,lmmaxvr
      do lm2=1,lmmaxvr
        do lm3=1,lmmaxvr
          zt1=gauntyry(lm2l(lm1),lm2l(lm2),lm2l(lm3),&
                   lm2m(lm1),lm2m(lm2),lm2m(lm3))
          if (abs(zt1).gt.1d-12) then
            do ir=1,nmtloc 
              zsummt=zsummt+zt1*rmtwt(ir)*dconjg(f1mt_(ir,lm1))*f2mt_(ir,lm2)*&
                f3mt_(ir,lm3)
            enddo
          endif
        enddo
      enddo
    enddo
  endif
! interstitial part
  do ir=1,ngrloc
    zsumir=zsumir+cfunir(ir+groffs)*dconjg(f1ir(ir,it))*f2ir(ir,it)*f3ir(ir,it)
  enddo
enddo
deallocate(f1mt_,f2mt_,f3mt_)
zt1=zsummt+zsumir*omega/dble(ngrtot)
call mpi_grid_reduce(zt1,all=.true.)
sic_int_zdz=zt1
end function

subroutine sic_mul_zd(alpha,f1mt,f1ir,f2mt,f2ir,f3mt,f3ir,tfmtuc)
use modmain
implicit none
! arguments
complex(8), intent(in) :: alpha
complex(8), intent(in) :: f1mt(lmmaxvr,nmtloc,sic_orbitals%ntr)
complex(8), intent(in) :: f1ir(ngrloc,sic_orbitals%ntr)
real(8), intent(in) :: f2mt(lmmaxvr,nmtloc,sic_orbitals%ntr)
real(8), intent(in) :: f2ir(ngrloc,sic_orbitals%ntr)
complex(8), intent(out) :: f3mt(lmmaxvr,nmtloc,sic_orbitals%ntr)
complex(8), intent(out) :: f3ir(ngrloc,sic_orbitals%ntr)
logical, intent(in) :: tfmtuc(sic_orbitals%ntr)
! local variables
complex(8), allocatable :: f1mt_(:,:),f3mt_(:,:)
real(8), allocatable :: f2mt_(:,:)
integer it,lm1,lm2,lm3,ir,lm
complex(8) zt1
complex(8), external :: gauntyry
allocate(f1mt_(nmtloc,lmmaxvr))
allocate(f2mt_(nmtloc,lmmaxvr))
allocate(f3mt_(nmtloc,lmmaxvr))
do it=1,sic_orbitals%ntr
  if (tfmtuc(it)) then
! muffin-tin part
    do lm=1,lmmaxvr
      f1mt_(:,lm)=f1mt(lm,:,it)
      f2mt_(:,lm)=f2mt(lm,:,it)
      f3mt_(:,lm)=zzero
    enddo
    do lm1=1,lmmaxvr
      do lm2=1,lmmaxvr
        do lm3=1,lmmaxvr
          zt1=gauntyry(lm2l(lm1),lm2l(lm2),lm2l(lm3),&
                   lm2m(lm1),lm2m(lm2),lm2m(lm3))
          if (abs(zt1).gt.1d-12) then
            do ir=1,nmtloc
              f3mt_(ir,lm3)=f3mt_(ir,lm3)+f1mt_(ir,lm1)*f2mt_(ir,lm2)*zt1
            enddo 
          endif
        enddo
      enddo
    enddo
    do lm=1,lmmaxvr
      f3mt(lm,:,it)=alpha*f3mt_(:,lm)
    enddo
  endif
  f3ir(:,it)=alpha*f1ir(:,it)*f2ir(:,it)
enddo
deallocate(f1mt_,f2mt_,f3mt_)
end subroutine

! <f|psi> = \int dr f^{*}(r) psi(r) =
!   = \sum_R \int_{Omega} dr f^{*}(r+R) psi(r+R) = 
!     \sum_R e^{ikR} \int_{Omega} dr f^{*}(r+R) psi(r)
complex(8) function sic_dot_lb(treduce,vpc,fmt1,fir1,tfmt1,fmt2,fir2)
use modmain
implicit none
logical, intent(in) :: treduce
real(8), intent(in) :: vpc(3)
complex(8), intent(in) :: fmt1(lmmaxvr,nmtloc,sic_orbitals%ntr)
complex(8), intent(in) :: fir1(ngrloc,sic_orbitals%ntr)
logical, intent(in) :: tfmt1(natmtot,sic_orbitals%ntr)
complex(8), intent(in) :: fmt2(lmmaxvr,nmtloc)
complex(8), intent(in) :: fir2(ngrloc)
complex(8) zdotmt,zdotir,zdot
complex(8), external :: zdotc
integer it,i,ias,ir
zdot=zzero
do it=1,sic_orbitals%ntr
  zdotmt=zzero
  zdotir=zzero
  do i=1,nmtloc
    ias=(mtoffs+i-1)/nrmtmax+1
    if (tfmt1(ias,it)) then
      zdotmt=zdotmt+rmtwt(i)*zdotc(lmmaxvr,fmt1(:,i,it),1,&
        fmt2(:,i),1)
    endif
  enddo
  do ir=1,ngrloc
    zdotir=zdotir+cfunir(ir+groffs)*dconjg(fir1(ir,it))*fir2(ir)
  enddo
  zdot=zdot+(zdotmt+zdotir*omega/dble(ngrtot))*&
    exp(zi*dot_product(vpc,sic_orbitals%vtc(:,it)))
enddo
if (treduce) call mpi_grid_reduce(zdot,all=.true.)
sic_dot_lb=zdot
return
end function

! this subroutine generates e^{igr} 
subroutine sic_genpw(vgpc,pwmt,pwir)
use modmain
implicit none
real(8), intent(in) :: vgpc(3)
complex(8), intent(out) :: pwmt(lmmaxvr,nmtloc)
complex(8), intent(out) :: pwir(ngrloc)
real(8) gpc
real(8) tpgp(2)
complex(8) ylmgp(lmmaxvr)
integer ias,is,ia,lm,ias1
real(8) jl(0:lmaxvr)
complex(8) zt1
integer ir,i
complex(8), allocatable :: pwmt_(:,:)

allocate(pwmt_(lmmaxvr,nmtloc))
pwmt_=zzero
! get spherical coordinates and length of G+q
call sphcrd(vgpc,gpc,tpgp)
! generate spherical harmonics for G+q
call genylm(lmaxvr,tpgp,ylmgp)

ias1=-1
do i=1,nmtloc
  ias=(mtoffs+i-1)/nrmtmax+1
  if (ias.ne.ias1) then
    is=ias2is(ias)
    ia=ias2ia(ias)
    zt1=fourpi*exp(zi*dot_product(vgpc,atposc(:,ia,is)))  
    ias1=ias
    ir=mod(mtoffs+i,nrmtmax)
  endif
  if (ir.le.nrmt(is)) then
! generate Bessel functions j_l(|G+q|x)
    call sbessel(lmaxvr,gpc*spr(ir,is),jl)
    do lm=1,lmmaxvr
      pwmt_(lm,i)=zt1*(zi**lm2l(lm))*jl(lm2l(lm))*dconjg(ylmgp(lm))
    enddo
  endif
  ir=ir+1
enddo
! convert to real spherical harmonics
!   remember that Y_{m}=\sum_{m'} dzsht_{m,m'} R_{m'}
call zgemm('T','N',lmmaxvr,nmtloc,lmmaxvr,zone,dzsht,lmmaxapw,pwmt_,lmmaxvr,&
  zzero,pwmt,lmmaxvr)
deallocate(pwmt_)
do ir=1,ngrloc
  pwir(ir)=exp(zi*dot_product(vgpc,vgrc(:,ir+groffs)))
enddo
return
end subroutine

subroutine sic_wavefmt(wffvmt,ist,fmt)
use modmain
implicit none
complex(8), intent(in) :: wffvmt(nstfv,lmmaxvr,nufrmax,natmtot)
integer, intent(in) :: ist
complex(8), intent(out) :: fmt(lmmaxvr,nmtloc)
integer ias,ias1,ir,is,ia,ic,i,lm,l,io
fmt=zzero

ias1=-1
do i=1,nmtloc
  ias=(mtoffs+i-1)/nrmtmax+1
  if (ias.ne.ias1) then
    is=ias2is(ias)
    ia=ias2ia(ias)
    ic=ias2ic(ias)
    ias1=ias
    ir=mod(mtoffs+i,nrmtmax)
  endif
  if (ir.le.nrmt(is)) then
    do lm=1,lmmaxvr
      l=lm2l(lm)
      do io=1,nufr(l,is)
        fmt(lm,i)=fmt(lm,i)+wffvmt(ist,lm,io,ias)*ufr(ir,l,io,ic)
      enddo
    enddo
  endif
  ir=ir+1
enddo

end subroutine



subroutine sic_gen_r(xmt,xir)
use modmain
implicit none
real(8), intent(out) :: xmt(lmmaxvr,nmtloc,sic_orbitals%ntr,4)
real(8), intent(out) :: xir(ngrloc,sic_orbitals%ntr,4)
integer ntp,it,itp,i,ias,ias1,ia,is,ir,lm
real(8), allocatable :: xmt_(:,:,:)
real(8), allocatable :: tp(:,:)
real(8), allocatable :: rlmb(:,:)
real(8), allocatable :: wtp(:)
real(8), allocatable :: spx(:,:)
real(8) rlm(lmmaxvr),vrc(3),t1

xmt=0.d0
xir=0.d0

ntp=266
allocate(tp(2,ntp))
allocate(rlmb(ntp,lmmaxvr))
allocate(wtp(ntp))
! Lebedev mesh
allocate(spx(3,ntp))
call leblaik(ntp,spx,wtp)                                                                   
do itp=1,ntp                   
  wtp(itp)=wtp(itp)*fourpi
  call sphcrd(spx(:,itp),t1,tp(:, itp))                                                     
enddo 
! generate spherical harmonics
do itp=1,ntp 
  call genrlm(lmaxvr,tp(1,itp),rlm)  
  do lm=1,lmmaxvr
    rlmb(itp,lm)=rlm(lm)*wtp(itp)
  enddo
enddo
deallocate(wtp,tp)

allocate(xmt_(ntp,nmtloc,4))
do it=1,sic_orbitals%ntr
  xmt_=0.d0
! compute r in muffin-tins
  ias1=-1
  do i=1,nmtloc
    ias=(mtoffs+i-1)/nrmtmax+1
    if (ias.ne.ias1) then
      is=ias2is(ias)
      ia=ias2ia(ias)
      ias1=ias
      ir=mod(mtoffs+i,nrmtmax)
    endif
    if (ir.le.nrmt(is)) then
      do itp=1,ntp
        vrc(:)=spx(:,itp)*spr(ir,is)+atposc(:,ia,is)+sic_orbitals%vtc(:,it)
        xmt_(itp,i,1)=vrc(1)
        xmt_(itp,i,2)=vrc(2)
        xmt_(itp,i,3)=vrc(3)
        xmt_(itp,i,4)=sum(vrc(:)**2)
      enddo
    endif
    ir=ir+1
  enddo !i
! expand in real spherical harmonics
  do i=1,4
    call dgemm('T','N',lmmaxvr,nmtloc,ntp,1.d0,rlmb,ntp,xmt_(1,1,i),ntp,0.d0,&
      xmt(1,1,it,i),lmmaxvr)
  enddo
! interstitial part
  do ir=1,ngrloc
    vrc(:)=vgrc(:,ir+groffs)+sic_orbitals%vtc(:,it)
    xir(ir,it,1)=vrc(1)
    xir(ir,it,2)=vrc(2)
    xir(ir,it,3)=vrc(3)
    xir(ir,it,4)=sum(vrc(:)**2)
  enddo
enddo !it
deallocate(rlmb,spx,xmt_)
end subroutine

end module
