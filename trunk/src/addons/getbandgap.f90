subroutine getbandgap(nkpt_,evalsv_,gap,ef)
use modmain
implicit none
integer, intent(in) :: nkpt_
real(8), intent(in) :: evalsv_(nstsv,nkpt_)
real(8), intent(out) :: gap
real(8), intent(out) :: ef
!
integer ist,i,j,nval
real(8) e(nstsv,2),e0,e1,e2(2)
! estimate the band gap
nval=nint(chgval)
gap=0.d0
if (spinpol.or.(.not.spinpol.and.mod(nval,2).eq.0)) then
  ist=nval
  if (.not.spinpol) ist=ist/2
  do j=1,nstsv
    e(j,1)=minval(evalsv_(j,:))
    e(j,2)=maxval(evalsv_(j,:))
  enddo
! sort bands using maximum band energies
  do i=1,nstsv-1
    do j=i+1,nstsv
      if (e(i,2).gt.e(j,2)) then
        e2(:)=e(i,:)
        e(i,:)=e(j,:)
        e(j,:)=e2(:)
      endif
    enddo
  enddo
  e0=e(ist,2)
  e1=e(ist+1,1)
  if (e1.gt.e0) then
    gap=e1-e0
    ef=e0+gap/2.d0
  endif  
endif
return
end subroutine
