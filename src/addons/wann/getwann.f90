subroutine getwann(ikloc)
use modmain
use mod_mpi_grid
implicit none
integer, intent(in) :: ikloc
integer ik_,nwann_,nspinor_,lmmaxvr_,nrfmax_,natmtot_,ngkmax_,recl
integer ik
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
inquire(iolength=recl)ik_,nwann_,nspinor_,lmmaxvr_,nrfmax_,natmtot_,ngkmax_, &
  wann_unkmt(:,:,:,:,:,ikloc),wann_unkit(:,:,:,ikloc),wann_c(:,:,ikloc)
open(70,file='WANN_UNK.OUT',action='READ',form='UNFORMATTED',access='DIRECT',recl=recl)
read(70,rec=ik)ik_,nwann_,nspinor_,lmmaxvr_,nrfmax_,natmtot_,ngkmax_, &
  wann_unkmt(:,:,:,:,:,ikloc),wann_unkit(:,:,:,ikloc),wann_c(:,:,ikloc)
close(70)
if (ik_.ne.ik.or.nwann_.ne.nwann.or.nspinor_.ne.nspinor.or. &
  lmmaxvr_.ne.lmmaxvr.or.nrfmax_.ne.nrfmax.or.natmtot_.ne.natmtot.or.ngkmax_.ne.ngkmax) then
  write(*,*)
  write(*,'("Error(getwann): wrong dimensions")')
  write(*,'("  mpi_grid_x : ",10I8)')mpi_grid_x
  write(*,'("  ik_ : ",I4," ik : ",I4)')ik_,ik
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
