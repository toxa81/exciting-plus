subroutine write_wann_ene
use modmain
implicit none
open(190,file="WANN_ENE.OUT",status="replace",form="unformatted")
write(190)wann_ene
close(190)
return
end