subroutine putwann(ik)
use modmain
implicit none
integer, intent(in) :: ik
integer recl
integer, external :: ikglob
inquire(iolength=recl)ik,nwann,nspinor,lmmaxvr,nrfmax,natmtot,ngkmax, &
  wann_unkmt(:,:,:,:,:,ik),wann_unkit(:,:,:,ik)
open(70,file='WANN_UNK.OUT',action='WRITE',form='UNFORMATTED', &
  access='DIRECT',recl=recl)
write(70,rec=ikglob(ik))ikglob(ik),nwann,nspinor,lmmaxvr,nrfmax,natmtot,ngkmax, &
  wann_unkmt(:,:,:,:,:,ik),wann_unkit(:,:,:,ik)
close(70)
return
end
