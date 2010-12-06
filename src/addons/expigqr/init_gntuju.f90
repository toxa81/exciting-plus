subroutine init_gntuju(iq,lmaxexp)
use modmain
use mod_wannier
use mod_expigqr
implicit none
integer, intent(in) :: iq
integer, intent(in) :: lmaxexp
call getmaxgnt(lmaxexp,ngntujumax)
if (allocated(ngntuju)) deallocate(ngntuju)
allocate(ngntuju(natmcls,ngvecme))
ngntuju=0
if (allocated(igntuju)) deallocate(igntuju)
allocate(igntuju(2,ngntujumax,natmcls,ngvecme))
igntuju=0
if (allocated(gntuju)) deallocate(gntuju)
allocate(gntuju(ngntujumax,natmcls,ngvecme))
gntuju=zzero
call gengntuju(iq,lmaxexp)
return
end