subroutine genwffvmt(lmax,lmmax,ngp,evecfv,apwalm,wffvmt)
use modmain
implicit none
! arguments
integer, intent(in) :: lmax
integer, intent(in) :: lmmax
integer, intent(in) :: ngp
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(out) :: wffvmt(nstfv,lmmax,nufrmax,natmtot)
! local variables
integer l,m,is,ia,ias,lm,i1,io,ilo,ig
integer ordl(0:lmax,natmtot)
! calculate first-variational coefficients
wffvmt=zzero
ordl=0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! apw coefficients
    do l=0,lmax
      do io=1,apword(l,is)
        ordl(l,ias)=ordl(l,ias)+1
        do m=-l,l
          lm=idxlm(l,m)
          call zgemv('T',ngp,nstfv,zone,evecfv,nmatmax, &
            apwalm(1,io,lm,ias),1,zzero,wffvmt(1,lm,ordl(l,ias),ias),1)
        enddo !m
      enddo !io
    enddo !l
! local orbital coefficients     
    do ilo=1,nlorb(is)
      l=lorbl(ilo,is)
      if (l.le.lmax) then
        ordl(l,ias)=ordl(l,ias)+1
        do m=-l,l
          lm=idxlm(l,m)
          i1=ngp+idxlo(lm,ilo,ias)
          wffvmt(:,lm,ordl(l,ias),ias)=evecfv(i1,:)
        enddo !m
      endif
    enddo !ilo    
  enddo
enddo
return
end