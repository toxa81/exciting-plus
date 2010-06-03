subroutine genwffvmt(lmax,lmmax,ngp,evecfv,apwalm,wffvmt)
use modmain
implicit none
! arguments
integer, intent(in) :: lmax
integer, intent(in) :: lmmax
integer, intent(in) :: ngp
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(out) :: wffvmt(lmmax,nufrmax,natmtot,nstfv)
! local variables
integer l,m,istfv,is,ia,ias,lm,i1,io,ilo,ig
integer ordl(0:lmax)
complex(8) zt1

wffvmt=zzero
! calculate first-variational coefficients
do istfv=1,nstfv
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      ordl=0
! apw coefficients
      do l=0,lmax
        do io=1,apword(l,is)
          ordl(l)=ordl(l)+1
          do m=-l,l
            lm=idxlm(l,m)
            zt1=dcmplx(0.d0,0.d0)
            do ig=1,ngp
              zt1=zt1+evecfv(ig,istfv)*apwalm(ig,io,lm,ias)
            enddo !ig
            wffvmt(lm,ordl(l),ias,istfv)=zt1
          enddo !m
        enddo !io
      enddo !l
! local orbital coefficients     
      do ilo=1,nlorb(is)
        l=lorbl(ilo,is)
        if (l.le.lmax) then
          ordl(l)=ordl(l)+1
          do m=-l,l
            lm=idxlm(l,m)
            i1=ngp+idxlo(lm,ilo,ias)
            wffvmt(lm,ordl(l),ias,istfv)=evecfv(i1,istfv)
          enddo !m
        endif
      enddo !ilo
    enddo !ia 
  enddo !is
enddo !istfv
return
end