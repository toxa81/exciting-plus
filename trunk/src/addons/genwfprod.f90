subroutine genwfprod(wfmt1,wfit1,wfmt2,wfit2,ngp1,ngp2,pwit,zprod)
use modmain
implicit none
! arguments
complex(8), intent(in) :: wfmt1(lmmaxvr,nrfmax,natmtot)
complex(8), intent(in) :: wfit1(ngkmax)
complex(8), intent(in) :: wfmt2(lmmaxvr,nrfmax,natmtot)
complex(8), intent(in) :: wfit2(ngkmax)
integer, intent(in) :: ngp1
integer, intent(in) :: ngp2
complex(8), intent(in) :: pwit(ngp1,ngp2)
complex(8), intent(out) :: zprod
! local variables
integer ig1,ig2,is,ia,ias,l,m,lm,io1,io2
zprod=zzero
! interstitial contribution
do ig1=1,ngp1
  do ig2=1,ngp2
    zprod=zprod+dconjg(wfit1(ig1))*pwit(ig1,ig2)*wfit2(ig2)
  enddo
enddo
! muffin-tin contribution
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do l=0,lmaxvr
      do io1=1,nrfmax
        do io2=1,nrfmax
          do m=-l,l
            lm=idxlm(l,m)
            zprod=zprod+dconjg(wfmt1(lm,io1,ias))*wfmt2(lm,io2,ias)*&
              urfprod(l,io1,io2,ias)
          enddo !m
        enddo
      enddo
    enddo !l
  enddo !ia
enddo !is
return
end
