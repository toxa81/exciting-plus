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
complex(8), intent(out) :: wfsvmt(lmmax,nufrmax,natmtot,nspinor,nstsv)
! local variables
integer ispn,j
complex(8), allocatable :: wffvmt(:,:,:,:)

allocate(wffvmt(nstfv,lmmax,nufrmax,natmtot))
call genwffvmt(lmax,lmmax,ngp,evecfv,apwalm,wffvmt)
! calculate second-variational coefficients
wfsvmt=zzero
do j=1,nstsv
  do ispn=1,nspinor
    call zgemv('T',nstfv,lmmax*nufrmax*natmtot,zone,wffvmt,nstfv,&
      evecsv(1+(ispn-1)*nstfv,j),1,zzero,wfsvmt(1,1,1,ispn,j),1)
  enddo
enddo
deallocate(wffvmt)
return
end
