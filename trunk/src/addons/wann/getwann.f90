subroutine getwann(ik)
use modmain
implicit none
integer, intent(in) :: ik
integer ik_,nwann_,nspinor_,lmmaxvr_,nrfmax_,natmtot_,ngkmax_,recl
integer, external :: ikglob
inquire(iolength=recl)ik_,nwann_,nspinor_,lmmaxvr_,nrfmax_,natmtot_,ngkmax_, &
  wann_unkmt(:,:,:,:,:,ik),wann_unkit(:,:,:,ik),wann_c(:,:,ik)
open(70,file='WANN_UNK.OUT',action='READ',form='UNFORMATTED',access='DIRECT',recl=recl)
read(70,rec=ikglob(ik))ik_,nwann_,nspinor_,lmmaxvr_,nrfmax_,natmtot_,ngkmax_, &
  wann_unkmt(:,:,:,:,:,ik),wann_unkit(:,:,:,ik),wann_c(:,:,ik)
close(70)
if (ik_.ne.ikglob(ik).or.nwann_.ne.nwann.or.nspinor_.ne.nspinor.or. &
  lmmaxvr_.ne.lmmaxvr.or.nrfmax_.ne.nrfmax.or.natmtot_.ne.natmtot.or.ngkmax_.ne.ngkmax) then
  write(*,*)
  write(*,'("Error(getwann): wrong dimensions")')
  write(*,'("  proc : ",I4)')iproc
  write(*,'("  ik_ : ",I4," ik : ",I4)')ik_,ikglob(ik)
  write(*,'("  nwann_ : ",I4," nwann : ",I4)')nwann_,nwann
  write(*,'("  nspinor_ : ",I4," nspinor : ",I4)')nspinor_,nspinor
  write(*,'("  lmmaxvr_ : ",I4," lmmaxvr : ",I4)')lmmaxvr_,lmmaxvr
  write(*,'("  nrfmax_ : ",I4," nrfmax : ",I4)')nrfmax_,nrfmax
  write(*,'("  natmtot_ : ",I4," natmtot : ",I4)')natmtot_,natmtot
  write(*,'("  ngkmax_ : ",I4," ngkmax : ",I4)')ngkmax_,ngkmax
  write(*,*)
  call pstop
endif
return
end
