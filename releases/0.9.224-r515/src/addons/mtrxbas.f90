subroutine mtrxbas(lmmax,tm,mtrx)
use modmain
integer, intent(in) :: lmmax
complex(8), intent(in) :: tm(lmmax,lmmax)
complex(8), intent(inout) :: mtrx(lmmax,lmmax,nspinor,nspinor)

integer lm1,lm2,lm3
complex(8), allocatable :: z1(:,:)

! transformation from Ylm(glob) to Rlm(loc): use rylm_lcs
! transformation from Rlm(loc) to Ylm(glob): use yrlm_lcs

allocate(z1(lmmax,lmmax))
do ispn=1,nspinor
  z1=zzero
  do lm1=1,lmmax
    do lm2=1,lmmax
      do lm3=1,lmmax
        z1(lm1,lm2)=z1(lm1,lm2)+dconjg(tm(lm1,lm3))*mtrx(lm3,lm2,ispn,ispn)
      enddo
    enddo
  enddo
  mtrx(1:lmmax,1:lmmax,ispn,ispn)=zzero
  do lm1=1,lmmax
    do lm2=1,lmmax
      do lm3=1,lmmax
        mtrx(lm1,lm2,ispn,ispn)=mtrx(lm1,lm2,ispn,ispn)+z1(lm1,lm3)*tm(lm2,lm3)
      enddo
    enddo
  enddo
enddo
deallocate(z1)
return
end
