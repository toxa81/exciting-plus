! this subroutine generates e^{ig(r+T)} 
subroutine genpw(vtl,vgpc,pwmt,pwir)
use modmain
implicit none
integer, intent(in) :: vtl(3)
real(8), intent(in) :: vgpc(3)
complex(8), intent(out) :: pwmt(lmmaxvr,nrmtmax,natmtot)
complex(8), intent(out) :: pwir(ngrtot)
real(8) :: vtrc(3) 
real(8) gpc
real(8) tpgp(2)
complex(8) ylmgp(lmmaxvr)
integer ias,is,ia,lm
real(8) jl(0:lmaxvr)
complex(8) zt1
integer ir,i1,i2,i3
real(8) v2(3),v3(3)

pwmt=zzero
pwir=zzero
vtrc(:)=vtl(1)*avec(:,1)+vtl(2)*avec(:,2)+vtl(3)*avec(:,3)
! get spherical coordinates and length of G+q
call sphcrd(vgpc,gpc,tpgp)
! generate spherical harmonics for G+q
call genylm(lmaxvr,tpgp,ylmgp)

do ias=1,natmtot
  is=ias2is(ias)
  ia=ias2ia(ias)
  zt1=fourpi*exp(zi*dot_product(vgpc,atposc(:,ia,is)+vtrc(:)))
  do ir=1,nrmt(is)
! generate Bessel functions j_l(|G+q|x)
    call sbessel(lmaxvr,gpc*spr(ir,is),jl)
    do lm=1,lmmaxvr
      pwmt(lm,ir,ias)=zt1*(zi**lm2l(lm))*jl(lm2l(lm))*dconjg(ylmgp(lm))
    enddo
  enddo
enddo
ir=0
do i3=0,ngrid(3)-1
  v2(3)=dble(i3)/dble(ngrid(3))
  do i2=0,ngrid(2)-1
    v2(2)=dble(i2)/dble(ngrid(2))
    do i1=0,ngrid(1)-1
      v2(1)=dble(i1)/dble(ngrid(1))
      call r3mv(avec,v2,v3)
      ir=ir+1
      pwir(ir)=exp(zi*dot_product(vgpc,v3(:)+vtrc(:))) 
    enddo
  enddo
enddo
return
end