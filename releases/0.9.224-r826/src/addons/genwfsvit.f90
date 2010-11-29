subroutine genwfsvit(ngp,evecfv,evecsv,wfsvit)
use modmain
implicit none
! arguments
integer, intent(in) :: ngp
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
complex(8), intent(out) :: wfsvit(ngkmax,nspinor,nstsv)
integer j,ispn,istfv
wfsvit=zzero
do j=1,nstsv
  do ispn=1,nspinor
    do istfv=1,nstfv
      wfsvit(1:ngp,ispn,j)=wfsvit(1:ngp,ispn,j)+&
        evecsv(istfv+(ispn-1)*nstfv,j)*evecfv(1:ngp,istfv)
    enddo
  enddo
enddo 
return
end
