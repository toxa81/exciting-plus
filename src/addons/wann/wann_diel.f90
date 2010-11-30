logical function wann_diel()
use modmain
use mod_wannier
implicit none
integer n
wann_diel=.true.
do n=1,nwantot
  if ((abs(wann_occ(n))*abs(wann_occ(n)-occmax)).gt.1d-8) wann_diel=.false.
enddo
return
end
