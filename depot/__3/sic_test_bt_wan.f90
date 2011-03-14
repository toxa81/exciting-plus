subroutine sic_test_bt_wan(n,wanlm)
use modmain
use mod_sic
implicit none
integer, intent(in) :: n
complex(8), intent(in) :: wanlm(s_ntp,s_nr,nspinor)
integer ir,itp
real(8) x(3),r0,t1
real(8), allocatable :: tp(:,:)
integer, parameter :: ntp=2000
integer, parameter :: nr=40
complex(8) zt1,wanval(nspinor)
character*100 fname

allocate(tp(2,ntp))
call sphcover(ntp,tp)

write(fname,'("wan_diff_.dat")')
open(209,file=trim(adjustl(fname)),form="formatted",status="replace")
do ir=1,nr
  write(fname,'("wan_diff_",I5.5,".dat")')ir
  open(210,file=trim(adjustl(fname)),form="formatted",status="replace")
  r0=sic_wan_cutoff*ir/dble(nr+1)
  write(210,'("# r : ",F12.6)')r0
  t1=0.d0
  do itp=1,ntp
    x(:)=r0*(/sin(tp(1,itp))*cos(tp(2,itp)),&
              sin(tp(1,itp))*sin(tp(2,itp)),&
              cos(tp(1,itp))/)
    zt1=s_func_val(x,wanlm(1,1,1))
    x(:)=x(:)-wanpos(:,n)
    call s_get_wanval(.true.,n,x,wanval)
    write(210,'(I8,G18.10," # ",2G18.10)')itp,abs(zt1-wanval(1)),dreal(zt1),dreal(wanval(1))
    t1=t1+abs(zt1-wanval(1))
  enddo
  write(210,*)
  write(210,'("# total diff : ",F12.6)')t1
  write(210,'("# average diff : ",F12.6)')t1/ntp
  close(210)
  write(209,'(2G18.10)')r0,t1
enddo
close(209)
return
end
