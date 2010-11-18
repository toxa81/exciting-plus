! this subroutine generates e^{igr} 
subroutine genpw(vgpc,pwmt,pwir)
use modmain
implicit none
real(8), intent(in) :: vgpc(3)
complex(8), intent(out) :: pwmt(lmmaxvr,nrmtmax,natmtot)
complex(8), intent(out) :: pwir(ngrtot)
real(8) gpc
real(8) tpgp(2)
complex(8) ylmgp(lmmaxvr)
integer ias,is,ia,lm
real(8) jl(0:lmaxvr)
complex(8) zt1
integer ir,i1,i2,i3,lm1,lm2
real(8) v2(3),v3(3)
complex(8) zt2(lmmaxvr)
complex(8), external :: zdotu

pwmt=zzero
pwir=zzero
! get spherical coordinates and length of G+q
call sphcrd(vgpc,gpc,tpgp)
! generate spherical harmonics for G+q
call genylm(lmaxvr,tpgp,ylmgp)

do ias=1,natmtot
  is=ias2is(ias)
  ia=ias2ia(ias)
  zt1=fourpi*exp(zi*dot_product(vgpc,atposc(:,ia,is)))
  do ir=1,nrmt(is)
! generate Bessel functions j_l(|G+q|x)
    call sbessel(lmaxvr,gpc*spr(ir,is),jl)
    do lm=1,lmmaxvr
      pwmt(lm,ir,ias)=zt1*(zi**lm2l(lm))*jl(lm2l(lm))*dconjg(ylmgp(lm))
    enddo
  enddo
enddo
! convert to real spherical harmonics
!   remember that Y_{m}=\sum_{m'} dzsht_{m,m'} R_{m'}
do ias=1,natmtot
  is=ias2is(ias)
  do ir=1,nrmt(is)
    zt2=zzero
    do lm1=1,lmmaxvr
      zt2(lm1)=zdotu(lmmaxvr,dzsht(1,lm1),1,pwmt(1,ir,ias),1)
    enddo
    pwmt(:,ir,ias)=zt2(:)
  enddo
enddo
do ir=1,ngrtot
  pwir(ir)=exp(zi*dot_product(vgpc,vgrc(:,ir)))
enddo
return
end