! this subroutine generates e^{igr} 
subroutine genpw(vgpc,pwmt,pwir,ylmtorlm)
use modmain
implicit none
real(8), intent(in) :: vgpc(3)
complex(8), intent(out) :: pwmt(lmmaxvr,nrmtmax,natmtot)
complex(8), intent(out) :: pwir(ngrtot)
complex(8), intent(in) :: ylmtorlm(lmmaxvr,lmmaxvr)
real(8) gpc
real(8) tpgp(2)
complex(8) ylmgp(lmmaxvr)
integer ias,is,ia,lm
real(8) jl(0:lmaxvr)
complex(8) zt1
integer ir,i1,i2,i3,lm1,lm2
real(8) v2(3),v3(3)
complex(8) zt2(lmmaxvr)

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
! Y_{m}=\sum_{m'} yrlm_{m,m'} R_{m'}
do ias=1,natmtot
  is=ias2is(ias)
  do ir=1,nrmt(is)
    zt2=zzero
    do lm1=1,lmmaxvr
      do lm2=1,lmmaxvr
        zt2(lm1)=zt2(lm1)+ylmtorlm(lm2,lm1)*pwmt(lm2,ir,ias)
      enddo
    enddo
    !pwmt(:,ir,ias)=zt2(:)
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
      pwir(ir)=exp(zi*dot_product(vgpc,v3(:))) 
    enddo
  enddo
enddo
return
end

subroutine genpwmt1(nvp,vpc,ias,ir,pwmt1)
use modmain
implicit none
integer, intent(in) :: nvp
real(8), intent(in) :: vpc(3,nvp)
integer, intent(in) :: ias
integer, intent(in) :: ir
complex(8), intent(out) :: pwmt1(nvp,lmmaxvr)
integer is,ia,ip,lm
real(8) pc,tpp(2),jl(0:lmaxvr)
complex(8) ylmp(lmmaxvr),zt1

is=ias2is(ias)
ia=ias2ia(ias)
do ip=1,nvp
  call sphcrd(vpc(1,ip),pc,tpp)
  call genylm(lmaxvr,tpp,ylmp)
  zt1=fourpi*exp(zi*dot_product(vpc(:,ip),atposc(:,ia,is)))
  call sbessel(lmaxvr,pc*spr(ir,is),jl)
  do lm=1,lmmaxvr
    pwmt1(ip,lm)=zt1*(zi**lm2l(lm))*jl(lm2l(lm))*dconjg(ylmp(lm))
  enddo
enddo
return
end

subroutine genpwir1(nvp,vpc,i3,i2,pwmt1)
use modmain
implicit none
integer, intent(in) :: nvp
real(8), intent(in) :: vpc(3,nvp)
integer, intent(in) :: i3
integer, intent(in) :: i2
complex(8), intent(out) :: pwmt1(nvp,ngrid(1))
real(8) v2(3),v3(3)
integer i1,ip
v2(3)=dble(i3)/dble(ngrid(3))
v2(2)=dble(i2)/dble(ngrid(2))
do i1=0,ngrid(1)-1
  v2(1)=dble(i1)/dble(ngrid(1))
  call r3mv(avec,v2,v3)
  do ip=1,nvp
    pwmt1(ip,i1+1)=exp(zi*dot_product(vpc(:,ip),v3(:))) 
  enddo
enddo
return
end


subroutine genpwmt2(ngpmax,ngp,igpig,ias,ir,pwmt2)
use modmain
implicit none
integer, intent(in) :: ngpmax
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngpmax)
integer, intent(in) :: ias
integer, intent(in) :: ir
complex(8), intent(out) :: pwmt2(ngpmax,lmmaxvr)
integer is,ia,igp,lm
real(8) gpc,tpgp(2),jl(0:lmaxvr)
complex(8) ylmgp(lmmaxvr),zt1

is=ias2is(ias)
ia=ias2ia(ias)
do igp=1,ngp
  call sphcrd(vgc(1,igpig(igp)),gpc,tpgp)
  call genylm(lmaxvr,tpgp,ylmgp)
  zt1=fourpi*exp(zi*dot_product(vgc(:,igpig(igp)),atposc(:,ia,is)))
  call sbessel(lmaxvr,gpc*spr(ir,is),jl)
  do lm=1,lmmaxvr
    pwmt2(igp,lm)=zt1*(zi**lm2l(lm))*jl(lm2l(lm))*dconjg(ylmgp(lm))
  enddo
enddo
return
end

subroutine genpwir2(ngpmax,ngp,igpig,i3,i2,pwmt2)
use modmain
implicit none
integer, intent(in) :: ngpmax
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngpmax)
integer, intent(in) :: i3
integer, intent(in) :: i2
complex(8), intent(out) :: pwmt2(ngpmax,ngrid(1))
real(8) v2(3),v3(3)
integer i1,igp
v2(3)=dble(i3)/dble(ngrid(3))
v2(2)=dble(i2)/dble(ngrid(2))
do i1=0,ngrid(1)-1
  v2(1)=dble(i1)/dble(ngrid(1))
  call r3mv(avec,v2,v3)
  do igp=1,ngp
    pwmt2(igp,i1+1)=exp(zi*dot_product(vgc(:,igpig(igp)),v3(:))) 
  enddo
enddo
return
end
