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
complex(8), intent(out) :: wfsvmt(lmmax,nrfmax,natmtot,nspinor,nstsv)
! local variables
integer ispn,istfv,j
complex(8), allocatable :: wffvmt(:,:,:,:)

allocate(wffvmt(lmmax,nrfmax,natmtot,nstfv))
call genwffvmt(lmax,lmmax,ngp,evecfv,apwalm,wffvmt)
! calculate second-variational coefficients
wfsvmt=zzero
do j=1,nstsv
  do ispn=1,nspinor
    do istfv=1,nstfv
      wfsvmt(:,:,:,ispn,j)=wfsvmt(:,:,:,ispn,j)+&
        evecsv(istfv+(ispn-1)*nstfv,j)*wffvmt(:,:,:,istfv)
    enddo
  enddo
enddo
deallocate(wffvmt)
return
end
