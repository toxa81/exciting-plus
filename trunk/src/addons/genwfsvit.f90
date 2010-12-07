subroutine genwfsvit(ngp,evecfv,evecsv,wfsvit)
use modmain
implicit none
! arguments
integer, intent(in) :: ngp
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
complex(8), intent(out) :: wfsvit(ngkmax,nspinor,nstsv)
integer j,ispn
wfsvit=zzero
do j=1,nstsv
  do ispn=1,nspinor
    call zgemv('N',ngp,nstfv,zone,evecfv,nmatmax,&
      evecsv(1+(ispn-1)*nstfv,j),1,zzero,wfsvit(1,ispn,j),1)
  enddo
enddo 
return
end
