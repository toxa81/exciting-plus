subroutine wfprodk(ngp,igpig,wfsvmt,wfsvit,wfnrmdev)
use modmain
implicit none
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
complex(8), intent(in) :: wfsvmt(lmmaxvr,nrfmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfsvit(ngkmax,nspinor,nstsv)
real(8), intent(out) :: wfnrmdev(nstsv*(nstsv+1)/2)

real(8) t1
integer ist1,ist2,j,ispn
complex(8) norm,zt1
complex(8), allocatable :: mit(:,:)

allocate(mit(ngp,ngp))

call genpwit(ngp,ngp,igpig,igpig,(/0,0,0/),mit)
j=0
do ist1=1,nstsv
  do ist2=ist1,nstsv
    j=j+1
    norm=zzero
    do ispn=1,nspinor
      call genwfprod(wfsvmt(1,1,1,ispn,ist1),wfsvit(1,ispn,ist1), &
        wfsvmt(1,1,1,ispn,ist2),wfsvit(1,ispn,ist2),ngp,ngp,mit,zt1)
      norm=norm+zt1
    enddo
    t1=0.d0
    if (ist1.eq.ist2) t1=1.d0
    wfnrmdev(j)=abs(norm-t1)
  enddo !ist1 
enddo !ist2

deallocate(mit)

return
end

