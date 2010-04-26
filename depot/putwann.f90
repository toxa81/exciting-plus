subroutine putwann(ikloc)
use modmain
use mod_mpi_grid
implicit none
integer, intent(in) :: ikloc
integer recl
integer :: ik
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
inquire(iolength=recl)ik,nwann,nspinor,lmmaxvr,nrfmax,natmtot,ngkmax, &
  wann_unkmt(:,:,:,:,:,ikloc),wann_unkit(:,:,:,ikloc),wann_c(:,:,ikloc)
open(70,file='WANN_UNK.OUT',action='WRITE',form='UNFORMATTED', &
  access='DIRECT',recl=recl)
write(70,rec=ik)ik,nwann,nspinor,lmmaxvr,nrfmax,natmtot,ngkmax, &
  wann_unkmt(:,:,:,:,:,ikloc),wann_unkit(:,:,:,ikloc),wann_c(:,:,ikloc)
close(70)
return
end
