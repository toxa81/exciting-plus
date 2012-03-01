!subroutine elk_init
!use modmain
!use mod_nrkp
!implicit none
!
!call readinput
!call init0
!call init1
!
!! read density and potentials from file
!call readstate
!! read Fermi energy from file
!call readfermi
!! generate the core wavefunctions and densities
!call gencore
!! find the new linearisation energies
!call linengy
!! generate the APW radial functions
!call genapwfr
!! generate the local-orbital radial functions
!call genlofr
!call getufr
!call genufrp
!wproc=.false.
!call genwfnr(-1,.false.,lmaxapw)
!
!return
!end
!
subroutine elk_m_init
use modmain
use mod_nrkp
use mod_madness
implicit none
integer ig,ik,ikloc
! allocate arrays for Madness-like WF generation (all MPI tasks
!  can resolve function at arbitrary point, no k-reduction is needed) 
if (allocated(m_ngknr)) deallocate(m_ngknr)
allocate(m_ngknr(nkptnr))
m_ngknr=0
if (allocated(m_igkignr)) deallocate(m_igkignr)
allocate(m_igkignr(ngkmax,nkptnr))
m_igkignr=0
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  m_ngknr(ik)=ngknr(ikloc)
  m_igkignr(:,ik)=igkignr(:,ikloc)
enddo
call mpi_grid_reduce(m_ngknr(1),nkptnr,dims=(/dim_k/),all=.true.)
call mpi_grid_reduce(m_igkignr(1,1),ngkmax*nkptnr,dims=(/dim_k/),all=.true.)
m_ngvec=0
do ik=1,nkptnr
  do ig=1,m_ngknr(ik)
    m_ngvec=max(m_ngvec,m_igkignr(ig,ik))
  enddo
enddo
if (allocated(m_wann_unkmt)) deallocate(m_wann_unkmt)
allocate(m_wann_unkmt(lmmaxapw,nufrmax,natmtot,nspinor,nkptnr))
if (allocated(m_wann_unkit)) deallocate(m_wann_unkit)
allocate(m_wann_unkit(ngkmax,nspinor,nkptnr))
return
end subroutine

subroutine elk_xc(dens,vx,vc,ex,ec)
use modmain
use modxcifc
implicit none
real(8), intent(in) :: dens(1,nspinor)
real(8), intent(out) :: vx(1,nspinor)
real(8), intent(out) :: vc(1,nspinor)
real(8), intent(out) :: ex(1)
real(8), intent(out) :: ec(1)

if (spinpol) then
  call xcifc(xctype,n=1,rhoup=dens(:,1),rhodn=dens(:,2),ex=ex(:),&
    &ec=ec(:),vxup=vx(:,1),vxdn=vx(:,2),vcup=vc(:,1),vcdn=vc(:,2))
else
  call xcifc(xctype,n=1,rho=dens(:,1),ex=ex(:),ec=ec(:),vx=vx(:,1),&
    &vc=vc(:,1))
endif
return
end subroutine

subroutine elk_load_wann_unk(n)
use modmain
use mod_nrkp
use mod_wannier
use mod_madness
implicit none
! arguments
integer, intent(in) :: n
! local variables
integer ikloc,ik,ias
!
m_wann_unkmt=zzero
m_wann_unkit=zzero
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  m_wann_unkmt(:,:,:,:,ik)=wann_unkmt(:,:,:,:,n,ikloc)
  m_wann_unkit(:,:,ik)=wann_unkit(:,:,n,ikloc)
enddo
call mpi_grid_reduce(m_wann_unkmt(1,1,1,1,1),&
  lmmaxapw*nufrmax*natmtot*nspinor*nkptnr,dims=(/dim_k/),all=.true.)
call mpi_grid_reduce(m_wann_unkit(1,1,1),&
  ngkmax*nspinor*nkptnr,dims=(/dim_k/),all=.true.)
m_wanpos(:)=wanpos(:,n)
return
end subroutine
!
!!subroutine elk_wann_unk_val(ik,n,ispn,r_cutoff,vrc,val)
!!use modmain
!!use mod_nrkp
!!use mod_wannier
!!! arguments
!!integer, intent(in) :: ik
!!integer, intent(in) :: n
!!integer, intent(in) :: ispn
!!real(8), intent(in) :: r_cutoff
!!real(8), intent(in) :: vrc(3)
!!real(8), intent(out) :: val(2)
!!! local variables
!!integer is,ia,ias,ir0,io,l,j,i,lm,ig
!!integer ntr(3)
!!real(8) vrc0(3),vtc(3),vr0(3),r0,tp(2),t1
!!real(8) ur(0:lmaxvr,nufrmax)
!!complex(8) zt1,zt2,zt3,ylm(lmmaxvr)
!!real(8) ya(nprad),c(nprad)
!!real(8), external :: polynom
!!logical, external :: vrinmt
!!
!!ias=wan_info(1,n)
!!vrc0(:)=vrc(:)-atposc(:,ias2ia(ias),ias2is(ias))
!!r0=sqrt(sum(vrc0(:)**2))
!!if (r0.gt.r_cutoff) then
!!  val(:)=0.d0
!!  return
!!endif
!!
!!zt1=zzero
!!zt2=zzero
!!if (vrinmt(vrc0,is,ia,ntr,vr0,ir0,r0)) then
!!  ias=idxas(ia,is)
!!  call sphcrd(vr0,t1,tp)
!!  call genylm(lmaxvr,tp,ylm)
!!  vtc(:)=ntr(1)*avec(:,1)+ntr(2)*avec(:,2)+ntr(3)*avec(:,3)
!!  ur=0.d0
!!  do l=0,lmaxvr
!!    do io=1,nufr(l,is)
!!      do j=1,nprad
!!        i=ir0+j-1
!!        ya(j)=ufr(i,l,io,ias2ic(ias))
!!      end do
!!      ur(l,io)=polynom(0,nprad,spr(ir0,is),ya,c,r0)
!!    enddo !io
!!  enddo !l
!!  zt3=exp(zi*dot_product(vkcnr(:,ik),vtc(:)))
!!  do lm=1,lmmaxvr
!!    l=lm2l(lm)
!!    do io=1,nufr(l,is)
!!      zt1=zt1+zt3*wann_unkmt1(lm,io,ias,ispn,n)*ur(l,io)*ylm(lm)
!!    enddo !io
!!  enddo !lm
!!else
!!  do ig=1,ngknr(ik)
!!    zt3=exp(zi*dot_product(vgkcnr(:,ig,ik),vrc(:)))/sqrt(omega)
!!    zt2=zt2+zt3*wann_unkit1(ig,ispn,n)
!!  enddo
!!endif
!!val(1)=dreal(zt1+zt2)
!!val(2)=dimag(zt1+zt2)
!!return
!!end
!
!
!
!subroutine elk_wan_val(n,ispn,r_cutoff,vrc,val)
!use modmain
!use mod_nrkp
!use mod_wannier
!implicit none
!integer, intent(in) :: n
!integer, intent(in) :: ispn
!real(8), intent(in) :: r_cutoff
!real(8), intent(in) :: vrc(3)
!real(8), intent(out) :: val(2)
!
!integer is,ia,ias,ir0,io,l,j,i,lm,ig
!integer ntr(3),ik
!real(8) vrc0(3),vtc(3),vr0(3),r0,tp(2),t1
!real(8) ur(0:lmaxvr,nufrmax)
!complex(8) zt1,zt2,zt3,ylm(lmmaxvr)
!real(8) ya(nprad),c(nprad)
!real(8), external :: polynom
!logical, external :: vrinmt
!
!ias=wan_info(1,n)
!vrc0(:)=vrc(:)-atposc(:,ias2ia(ias),ias2is(ias))
!r0=sqrt(sum(vrc0(:)**2))
!if (r0.gt.r_cutoff) then
!  val(:)=0.d0
!  return
!endif
!
!zt1=zzero
!zt2=zzero
!if (vrinmt(vrc0,is,ia,ntr,vr0,ir0,r0)) then
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
!  do ik=1,nkptnr
!    zt3=exp(zi*dot_product(vkcnr(:,ik),vtc(:)))
!    do lm=1,lmmaxvr
!      l=lm2l(lm)
!      do io=1,nufr(l,is)
!        zt1=zt1+zt3*wann_unkmt(lm,io,ias,ispn,n,ik)*ur(l,io)*ylm(lm)
!      enddo !io
!    enddo !lm
!  enddo !ik
!else
!  do ik=1,nkptnr
!    do ig=1,ngknr(ik)
!      zt3=exp(zi*dot_product(vgkcnr(:,ig,ik),vrc(:)))/sqrt(omega)
!      zt2=zt2+zt3*wann_unkit(ig,ispn,n,ik)
!    enddo
!  enddo
!endif
!val(1)=dreal(zt1+zt2)/nkptnr
!val(2)=dimag(zt1+zt2)/nkptnr
!return
!end
!
!
subroutine elk_wan_rho(x,rcutoff,wrho)
use modmain
use mod_nrkp
use mod_wannier
use mod_madness
use mod_sic
implicit none
real(8), intent(inout) :: x(3)
real(8), intent(in) :: rcutoff
real(8), intent(out) :: wrho
complex(8) wanval(nspinor)
integer ispn
!
wrho=0.d0
call s_get_wanval(x,wanval,rcutoff=rcutoff)
do ispn=1,nspinor
  wrho=wrho+abs(wanval(ispn))**2
enddo
return
end subroutine

subroutine sic_wan_rho(j,x,rcutoff,wrho)
use modmain
use mod_nrkp
use mod_wannier
use mod_madness
use mod_sic
implicit none
integer, intent(in) :: j
real(8), intent(inout) :: x(3)
real(8), intent(in) :: rcutoff
real(8), intent(out) :: wrho
complex(8) wanval(nspinor)
integer ispn
!
wrho=0.d0
call s_spinor_func_val(x,s_wlm(1,1,1,j),wanval,rcutoff=rcutoff)
do ispn=1,nspinor
  wrho=wrho+abs(wanval(ispn))**2
enddo
return
end subroutine


