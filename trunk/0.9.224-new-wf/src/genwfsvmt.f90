subroutine genwfsvmt(lmax,lmmax,ngp,evecfv,evecsv,apwalm,wfsvmt)
use modmain
implicit none
! arguments
integer, intent(in) :: lmax
integer, intent(in) :: lmmax
integer, intent(in) :: ngp
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(out) :: wfsvmt(lmmax,nrfmax,natmtot,nstsv,nspinor)
! local variables
integer j,l,m,ispn,istfv,is,ia,ias,lm,ig,i1,io,ilo,iwf
integer ordl(0:lmax)
complex(8), allocatable :: wffvmt(:,:,:,:)
complex(8) zt1

allocate(wffvmt(nstfv,nrfmax,lmmax,natmtot))
call genwffvmt(lmax,lmmax,ngp,evecfv,apwalm,wffvmt)
! calculate second-variational coefficients
wfsvmt=dcmplx(0.d0,0.d0)
do iwf=1,nstsv
  do ias=1,natmtot
    do io=1,nrfmax
      do lm=1,lmmax
        do ispn=1,nspinor
          zt1=dcmplx(0.d0,0.d0)
          do istfv=1,nstfv
            zt1=zt1+evecsv(istfv+(ispn-1)*nstfv,iwf)*wffvmt(istfv,io,lm,ias)
          enddo !istfv
          wfsvmt(lm,io,ias,iwf,ispn)=zt1
        enddo !ispn
      enddo !lm
    enddo !io
  enddo !ias
enddo !iwf
deallocate(wffvmt)
return
end

subroutine genwffvmt(lmax,lmmax,ngp,evecfv,apwalm,wffvmt)
use modmain
implicit none
! arguments
integer, intent(in) :: lmax
integer, intent(in) :: lmmax
integer, intent(in) :: ngp
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(out) :: wffvmt(nstfv,nrfmax,lmmax,natmtot)
! local variables
integer l,m,istfv,is,ia,ias,lm,i1,io,ilo,ig
integer ordl(0:lmax)
complex(8) zt1

wffvmt=dcmplx(0.d0,0.d0)
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
            wffvmt(istfv,ordl(l),lm,ias)=zt1
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
            wffvmt(istfv,ordl(l),lm,ias)=evecfv(i1,istfv)
          enddo !m
        endif
      enddo !ilo
    enddo !ia 
  enddo !is
enddo !istfv

return
end

