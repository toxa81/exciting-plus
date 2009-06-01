subroutine genwann_c(ik,e,wfsvmt,wann_c_)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(in) :: e(nstsv)
complex(8), intent(in) :: wfsvmt(lmmaxvr,nrfmax,natmtot,nspinor,nstsv)
complex(8), intent(out) :: wann_c_(nwann,nstsv)
! local variables
complex(8), allocatable :: prjao(:,:)
complex(8), allocatable :: s(:,:),sdiag(:)
integer ispn,j,n,m1,m2,io1,io2,ias,lm1,lm,ierr,l,itype
logical l1

! compute <\psi|g_n>
allocate(prjao(nwann,nstsv))
prjao=zzero
do n=1,nwann
  ias=iwann(1,n)
  lm=iwann(2,n)
  l=lm2l(lm)
  ispn=iwann(3,n)
  itype=iwann(4,n)
  do j=1,nstsv
    l1=.false.
    if (wann_use_eint) then
      if (e(j).ge.wann_eint(1,itype).and.e(j).le.wann_eint(2,itype)) l1=.true.
    else
      if (j.ge.wann_nint(1,itype).and.j.le.wann_nint(2,itype)) l1=.true.
    endif
    if (l1) then
      do m1=-l,l
        lm1=idxlm(l,m1)
        io2=2
        do io1=1,nrfmax
          prjao(n,j)=prjao(n,j)+dconjg(wfsvmt(lm1,io1,ias,ispn,j)) * &
            urfprod(l,io1,io2,ias)*rylm_lcs(lm,lm1,ias)
        enddo !io1
      enddo !m
      if (abs(prjao(n,j)).lt.1d-2) prjao(n,j)=zzero
    endif
  enddo !j
enddo !n
! compute ovelap matrix
allocate(s(nwann,nwann))
allocate(sdiag(nwann))
s=zzero
do m1=1,nwann
  do m2=1,nwann
    do j=1,nstsv
      s(m1,m2)=s(m1,m2)+prjao(m1,j)*dconjg(prjao(m2,j))
    enddo
  enddo
  sdiag(m1)=s(m1,m1)
enddo
! compute S^{-1/2}
call isqrtzhe(nwann,s,ierr)
if (ierr.ne.0) then
  write(*,*)
  write(*,'("Warning(genwann_c): failed to calculate S^{-1/2}")')
  write(*,'("  k-point : ",I4)')ik
  write(*,'("  iteration : ",I4)')iscl
  write(*,'("  number of linear dependent WFs : ",I4)')ierr
  write(*,'("  diagonal elements of overlap matrix : ")')
  write(*,'(6X,5G18.10)')abs(sdiag)
  write(*,'("Non-orthogonal WFs will be used")')
  write(*,*)
endif
! compute Wannier function expansion coefficients
wann_c_=zzero
if (ierr.eq.0) then
  do m1=1,nwann
    do m2=1,nwann
      wann_c_(m1,:)=wann_c_(m1,:)+prjao(m2,:)*dconjg(s(m2,m1))
    enddo
  enddo
else
  wann_c_=prjao
endif
deallocate(s,sdiag)
deallocate(prjao)

return
end
