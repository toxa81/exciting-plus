subroutine unimtrxt(lmmax,tm,mtrx)
use modmain
implicit none
integer, intent(in) :: lmmax
complex(8), intent(in) :: tm(lmmax,lmmax)
complex(8), intent(inout) :: mtrx(lmmax,lmmax)
integer lm1,lm2,lm3
complex(8), allocatable :: z1(:,:)
! transformation from Ylm(glob) to Rlm(loc): use rylm_lps
! transformation from Rlm(loc) to Ylm(glob): use yrlm_lps
! TODO: convert to blas calls
allocate(z1(lmmax,lmmax))
z1=zzero
do lm1=1,lmmax
  do lm2=1,lmmax
    do lm3=1,lmmax
      z1(lm1,lm2)=z1(lm1,lm2)+dconjg(tm(lm1,lm3))*mtrx(lm3,lm2)
    enddo
  enddo
enddo
mtrx(1:lmmax,1:lmmax)=zzero
do lm1=1,lmmax
  do lm2=1,lmmax
    do lm3=1,lmmax
      mtrx(lm1,lm2)=mtrx(lm1,lm2)+z1(lm1,lm3)*tm(lm2,lm3)
    enddo
  enddo
enddo
deallocate(z1)
return
end
