subroutine print_info
use modmain
implicit none

write(*,'(80("="))')
write(*,'(" task : ",I4)')task
write(*,'(80("="))')
write(*,'(" ngvec",T50,":",I6)')ngvec
write(*,'(" gkmax",T50,":",G18.10)')gkmax
write(*,'(" ngkmax",T50,":",I6)')ngkmax
write(*,'(" nlotot",T50,":",I6)')nlotot
write(*,'(" nmatmax",T50,":",I6)')nmatmax

return
end
