subroutine genprjao(ias,lm,ispn,i,wfsvmt,prjao)
use modmain
implicit none
integer, intent(in) :: ias
integer, intent(in) :: lm
integer, intent(in) :: ispn
integer, intent(in) :: i
complex(8), intent(in) :: wfsvmt(lmmaxvr,nufrmax,natmtot,nspinor,nstsv)
complex(8), intent(out) :: prjao

integer l,m1,lm1,io1,ir,is,ic
real(8) fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)
! compute <psi_{ik}|phi_n>, where n={ias,lm,ispn} 
! |psi> is a spinor Bloch-function 
! |phi> is a valence local orbital
! 
l=lm2l(lm)
is=ias2is(ias)
ic=ias2ic(ias)
prjao=zzero
do m1=-l,l
  lm1=idxlm(l,m1)
  do io1=1,nufr(l,is)
! project to local orbital    
    if (wannier_prjao.eq.0) then
      prjao=prjao+dconjg(wfsvmt(lm1,io1,ias,ispn,i))*&
        ufrp(l,io1,apword(l,is)+1,ic)*rylm_lps(lm,lm1,ias)
    endif
! project to f(x)=(1+cos(Pi*x/R))
    if (wannier_prjao.eq.1) then
      do ir=1,nrmt(is)
        fr(ir)=ufr(ir,l,io1,ic)*(1+cos(pi*spr(ir,is)/rmt(is)))*(spr(ir,is)**2)                                                        
      enddo
      call fderiv(-1,nrmt(is),spr(1,is),fr,gr,cf)
      prjao=prjao+dconjg(wfsvmt(lm1,io1,ias,ispn,i))*gr(nrmt(is))*rylm_lps(lm,lm1,ias)
    endif
  enddo !io1
enddo !m
return
end
