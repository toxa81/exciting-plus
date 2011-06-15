subroutine genwfsvocc(lmax,lmmax,ngp,nstocc,istocc,evecfv,evecsv,&
  apwalm,wfsvmt,wfsvit)
use modmain
implicit none
! arguments
integer, intent(in) :: lmax
integer, intent(in) :: lmmax
integer, intent(in) :: ngp
integer, intent(in) :: nstocc
integer, intent(in) :: istocc(nstsv)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(out) :: wfsvmt(lmmax,nufrmax,natmtot,nspinor,nstocc)
complex(8), intent(out) :: wfsvit(ngkmax,nspinor,nstocc) 
! local variables
integer ispn,nrow,ist,jst
complex(8), allocatable :: wffvmt(:,:,:,:)
!
allocate(wffvmt(lmmax,nufrmax,natmtot,nstfv))
call genwffvmt(lmax,lmmax,ngp,evecfv,apwalm,wffvmt)
nrow=lmmax*nufrmax*natmtot
wfsvmt=zzero
wfsvit=zzero
do jst=1,nstocc
  ist=istocc(jst)
  do ispn=1,nspinor
    call zgemv('N',nrow,nstfv,zone,wffvmt,nrow,&
      evecsv((ispn-1)*nstfv+1,ist),1,zzero,wfsvmt(1,1,1,ispn,jst),1)
    call zgemv('N',ngp,nstfv,zone,evecfv,nmatmax,&
      evecsv((ispn-1)*nstfv+1,ist),1,zzero,wfsvit(1,ispn,jst),1)
  enddo
enddo
deallocate(wffvmt)
return
end subroutine
