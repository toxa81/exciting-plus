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
