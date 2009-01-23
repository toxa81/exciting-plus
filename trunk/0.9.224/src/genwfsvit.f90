subroutine genwfsvit(ngp,evecfv,evecsv,wfsvit)
use modmain
implicit none
! arguments
integer, intent(in) :: ngp
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
complex(8), intent(out) :: wfsvit(ngkmax,nstsv,nspinor)

integer iwf,ig,ispn,istfv
complex(8) zt1

wfsvit=dcmplx(0.d0,0.d0)
do iwf=1,nstsv
  do ig=1,ngp
    do ispn=1,nspinor
      zt1=dcmplx(0.d0,0.d0)
      do istfv=1,nstfv
        zt1=zt1+evecsv(istfv+(ispn-1)*nstfv,iwf)*evecfv(ig,istfv)
      enddo
      wfsvit(ig,iwf,ispn)=zt1
    enddo 
  enddo !ig
enddo !iwf

return
end
