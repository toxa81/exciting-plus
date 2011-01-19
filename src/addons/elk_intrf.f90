subroutine elk_init
use modmain
use mod_nrkp
implicit none

call readinput
call init0
call init1

! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! generate the core wavefunctions and densities
call gencore
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
call getufr
call genufrp
wproc=.false.
call genwfnr(-1,.false.)

return
end

subroutine elk_xc(dens,vxc,exc)
use modmain
use modxcifc
implicit none
real(8), intent(in) :: dens(nspinor)
real(8), intent(out) :: vxc(nspinor)
real(8), intent(out) :: exc
! local variables
real(8) ex(1),ec(1),vx1(1),vc1(1),vx2(1),vc2(1)

if (spinpol) then
  call xcifc(xctype,n=1,rhoup=dens(1),rhodn=dens(2),ex=ex,ec=ec,&
    vxup=vx1,vxdn=vx2,vcup=vc1,vcdn=vc2)
  vxc(1)=vx1(1)+vc1(1)
  vxc(2)=vx2(1)+vc2(1)
else
  call xcifc(xctype,n=1,rho=dens(1),ex=ex,ec=ec,vx=vx1,vc=vc1)
  vxc(1)=vx1(1)+vc1(1)
endif
exc=ex(1)+ec(1)
return
end subroutine

subroutine elk_load_wann_unk(ik)
use modmain
use mod_nrkp
use mod_wannier
implicit none
!
integer, intent(in) :: ik
!
integer ikloc,h
ikloc=mpi_grid_map(nkptnr,dim_k,glob=ik,x=h)
if (mpi_grid_dim_pos(dim_k).eq.h) then
  wann_unkmt1(:,:,:,:,:)=wann_unkmt(:,:,:,:,:,ikloc)
  wann_unkit1(:,:,:)=wann_unkit(:,:,:,ikloc)
endif
call mpi_grid_bcast(wann_unkmt1(1,1,1,1,1),&
  lmmaxvr*nufrmax*natmtot*nspinor*nwantot,root=(/h/))
call mpi_grid_bcast(wann_unkit1(1,1,1),&
  ngkmax*nspinor*nwantot,root=(/h/))
return
end subroutine

subroutine elk_wann_unk_val(ik,n,ispn,r_cutoff,vrc,val)
use modmain
use mod_nrkp
use mod_wannier
! arguments
integer, intent(in) :: ik
integer, intent(in) :: n
integer, intent(in) :: ispn
real(8), intent(in) :: r_cutoff
real(8), intent(in) :: vrc(3)
real(8), intent(out) :: val(2)
! local variables
integer is,ia,ias,ir0,io,l,j,i,lm,ig
integer ntr(3)
real(8) vrc0(3),vtc(3),vr0(3),r0,tp(2),t1
real(8) ur(0:lmaxvr,nufrmax)
complex(8) zt1,zt2,zt3,ylm(lmmaxvr)
real(8) ya(nprad),c(nprad)
real(8), external :: polynom
logical, external :: vrinmt

ias=wan_info(1,n)
vrc0(:)=vrc(:)-atposc(:,ias2ia(ias),ias2is(ias))
r0=sqrt(sum(vrc0(:)**2))
if (r0.gt.r_cutoff) then
  val(:)=0.d0
  return
endif

zt1=zzero
zt2=zzero
if (vrinmt(vrc0,is,ia,ntr,vr0,ir0,r0)) then
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
  zt3=exp(zi*dot_product(vkcnr(:,ik),vtc(:)))
  do lm=1,lmmaxvr
    l=lm2l(lm)
    do io=1,nufr(l,is)
      zt1=zt1+zt3*wann_unkmt1(lm,io,ias,ispn,n)*ur(l,io)*ylm(lm)
    enddo !io
  enddo !lm
else
  do ig=1,ngknr(ik)
    zt3=exp(zi*dot_product(vgkcnr(:,ig,ik),vrc(:)))/sqrt(omega)
    zt2=zt2+zt3*wann_unkit1(ig,ispn,n)
  enddo
endif
val(1)=dreal(zt1+zt2)
val(2)=dimag(zt1+zt2)
return
end



subroutine elk_wan_val(n,ispn,r_cutoff,vrc,val)
use modmain
use mod_nrkp
use mod_wannier
implicit none
integer, intent(in) :: n
integer, intent(in) :: ispn
real(8), intent(in) :: r_cutoff
real(8), intent(in) :: vrc(3)
real(8), intent(out) :: val(2)

integer is,ia,ias,ir0,io,l,j,i,lm,ig
integer ntr(3),ik
real(8) vrc0(3),vtc(3),vr0(3),r0,tp(2),t1
real(8) ur(0:lmaxvr,nufrmax)
complex(8) zt1,zt2,zt3,ylm(lmmaxvr)
real(8) ya(nprad),c(nprad)
real(8), external :: polynom
logical, external :: vrinmt

ias=wan_info(1,n)
vrc0(:)=vrc(:)-atposc(:,ias2ia(ias),ias2is(ias))
r0=sqrt(sum(vrc0(:)**2))
if (r0.gt.r_cutoff) then
  val(:)=0.d0
  return
endif

zt1=zzero
zt2=zzero
if (vrinmt(vrc0,is,ia,ntr,vr0,ir0,r0)) then
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
  do ik=1,nkptnr
    zt3=exp(zi*dot_product(vkcnr(:,ik),vtc(:)))
    do lm=1,lmmaxvr
      l=lm2l(lm)
      do io=1,nufr(l,is)
        zt1=zt1+zt3*wann_unkmt(lm,io,ias,ispn,n,ik)*ur(l,io)*ylm(lm)
      enddo !io
    enddo !lm
  enddo !ik
else
  do ik=1,nkptnr
    do ig=1,ngknr(ik)
      zt3=exp(zi*dot_product(vgkcnr(:,ig,ik),vrc(:)))/sqrt(omega)
      zt2=zt2+zt3*wann_unkit(ig,ispn,n,ik)
    enddo
  enddo
endif
val(1)=dreal(zt1+zt2)/nkptnr
val(2)=dimag(zt1+zt2)/nkptnr
return
end


subroutine elk_wan_rho(n,r_cutoff,vrc,wrho)
use modmain
use mod_nrkp
use mod_wannier
implicit none
integer, intent(in) :: n
real(8), intent(in) :: r_cutoff
real(8), intent(in) :: vrc(3)
real(8), intent(out) :: wrho
!
integer is,ia,ias,ir0,io,l,j,i,lm,ig,ispn
integer ntr(3),ik
real(8) vrc0(3),vtc(3),vr0(3),r0,tp(2),t1
real(8) ur(0:lmaxvr,nufrmax)
complex(8) wanval(nspinor),zt1,zt2,zt3,ylm(lmmaxvr)
real(8) ya(nprad),c(nprad)
real(8), external :: polynom
logical, external :: vrinmt
complex(8) expigr(ngvec)
complex(8) zm1(lmmaxvr,nufrmax,nspinor)

wrho=0.d0

ias=wan_info(1,n)
vrc0(:)=vrc(:)-atposc(:,ias2ia(ias),ias2is(ias))
r0=sqrt(sum(vrc0(:)**2))
if (r0.gt.r_cutoff) return

wanval=zzero
if (vrinmt(vrc0,is,ia,ntr,vr0,ir0,r0)) then
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
!  do ik=1,nkptnr
!    zt1=exp(zi*dot_product(vkcnr(:,ik),vtc(:)))
!    do ispn=1,nspinor
!      zt2=zzero
!      do lm=1,lmmaxvr
!        l=lm2l(lm)
!        do io=1,nufr(l,is)
!          zt2=zt2+wann_unkmt(lm,io,ias,ispn,n,ik)*ur(l,io)*ylm(lm)
!        enddo !io
!      enddo !lm
!      wanval(ispn)=wanval(ispn)+zt2*zt1
!    enddo !ispn
!  enddo !ik
  zm1=zzero
  do ik=1,nkptnr
    zt1=exp(zi*dot_product(vkcnr(:,ik),vtc(:)))
    do ispn=1,nspinor
      zm1(:,:,ispn)=zm1(:,:,ispn)+zt1*wann_unkmt(:,:,ias,ispn,n,ik)
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
  do ig=1,ngvec
    expigr(ig)=exp(zi*dot_product(vrc(:),vgc(:,ig)))  
  enddo
  do ik=1,nkptnr
    zt1=exp(zi*dot_product(vrc(:),vkcnr(:,ik)))/sqrt(omega)
    do ig=1,ngknr(ik)
      zt2=expigr(igkignr(ig,ik))
      do ispn=1,nspinor
        wanval(ispn)=wanval(ispn)+zt1*zt2*wann_unkit(ig,ispn,n,ik)
      enddo
    enddo
  enddo
endif
do ispn=1,nspinor
  wrho=wrho+abs(wanval(ispn)/nkptnr)**2
enddo
return
end
