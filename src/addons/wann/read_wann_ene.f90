subroutine read_wann_ene
use modmain
implicit none
open(190,file="WANN_ENE.OUT",status="old",form="unformatted")
read(190)wann_ene
close(190)
return
end