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
integer ispn,j,nrow
complex(8), allocatable :: wffvmt(:,:,:,:)
complex(8), allocatable :: wfsvmt_(:,:,:,:)
!
allocate(wffvmt(lmmax,nufrmax,natmtot,nstfv))
call genwffvmt(lmax,lmmax,ngp,evecfv,apwalm,wffvmt)
allocate(wfsvmt_(lmmax,nufrmax,natmtot,nstsv))
nrow=lmmax*nufrmax*natmtot
! calculate second-variational coefficients
wfsvmt=zzero
do ispn=1,nspinor
  call zgemm('N','N',nrow,nstsv,nstfv,zone,wffvmt,nrow,&
    evecsv((ispn-1)*nstfv+1,1),nstsv,zzero,wfsvmt_,nrow)
  do j=1,nstsv
    wfsvmt(:,:,:,ispn,j)=wfsvmt_(:,:,:,j) 
  enddo
enddo
deallocate(wffvmt,wfsvmt_)
return
end
