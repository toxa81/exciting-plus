subroutine genpwit(ngp1,ngp2,igpig1,igpig2,vgl3,pwit)
use modmain
implicit none
! arguments
integer, intent(in) :: ngp1
integer, intent(in) :: ngp2
integer, intent(in) :: igpig1(ngkmax)
integer, intent(in) :: igpig2(ngkmax)
integer, intent(in) :: vgl3(3)
complex(8), intent(out) :: pwit(ngp1,ngp2)
! local variables
integer vecgl(3),ig1,ig2,is,ia,ias
real(8) vecgc(3),vtmp(3),lengc,tpgc(2)
logical l0
complex(8) sfacgc(natmtot)

pwit=zzero
do ig1=1,ngp1
  do ig2=1,ngp2
! -G1+G2+G3 vector in lattice coordinates
    vecgl(:)=-ivg(:,igpig1(ig1))+ivg(:,igpig2(ig2))+vgl3(:)
    l0=.false.
    if (sum(abs(vecgl)).eq.0) then
      pwit(ig1,ig2)=zone
      l0=.true.
    endif
    vtmp(:)=1.d0*vecgl(:)
    call r3mv(bvec,vtmp,vecgc)
    call sphcrd(vecgc,lengc,tpgc)
    call gensfacgp(1,vecgc,1,sfacgc)
    do is=1,nspecies
      do ia=1,natoms(is)
        ias=idxas(ia,is)
        if (l0) then
          pwit(ig1,ig2)=pwit(ig1,ig2)-(fourpi/omega)*sfacgc(ias)*&
            (rmt(is)**3)/3.d0
        else
          pwit(ig1,ig2)=pwit(ig1,ig2)-(fourpi/omega)*sfacgc(ias)*&
            (-lengc*rmt(is)*cos(lengc*rmt(is))+sin(lengc*rmt(is)))/lengc**3
        endif
      enddo !ia
    enddo !is
  enddo
enddo

return
end
