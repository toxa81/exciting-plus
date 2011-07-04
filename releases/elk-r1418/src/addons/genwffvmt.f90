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
integer l,m,is,ia,ias,lm,i1,io,ilo,nrow
integer ordl(0:lmax,natmtot)
complex(8), allocatable :: apwalm_(:,:,:,:)
complex(8), allocatable :: wffvmt_(:,:,:,:)
! change order of indexes
allocate(apwalm_(ngkmax,lmmax,apwordmax,natmtot))
do ias=1,natmtot
  do io=1,apwordmax
    do lm=1,lmmax
      apwalm_(:,lm,io,ias)=apwalm(:,io,lm,ias)
    enddo
  enddo
enddo
! precompute \sum_{G} apwalm_(G,lm,io,ias)*evecfv(G,ist)
allocate(wffvmt_(lmmax,apwordmax,natmtot,nstfv))
nrow=lmmax*apwordmax*natmtot
call zgemm('T','N',nrow,nstfv,ngp,zone,apwalm_,ngkmax,evecfv,nmatmax,zzero,&
  wffvmt_,nrow)
! get first-variational coefficients
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
          wffvmt(lm,ordl(l,ias),ias,:)=wffvmt_(lm,ordl(l,ias),ias,:)
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
          wffvmt(lm,ordl(l,ias),ias,:)=evecfv(i1,:)
        enddo !m
      endif
    enddo !ilo    
  enddo
enddo
deallocate(apwalm_,wffvmt_)
return
end