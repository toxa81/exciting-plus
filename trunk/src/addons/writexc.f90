subroutine writexc
use modmain
use modxcifc
implicit none
integer i
real(8) d1(1),d2(1),d3(1),d4(1),d5(1)
call init0
call init1

do i=1,100
  d1=i/100.d0
  call xcifc(xctype,n=1,rho=d1,ex=d2,ec=d3,vx=d4,vc=d5)
  write(100,*)d1(1),d4(1)+d5(1)
enddo
end

!call xcifc(xctype,n=ngrtot,rhoup=rfir(:,1),rhodn=rfir(:,2),ex=exir, &
!     ec=ecir,vxup=vx(:,1),vxdn=vx(:,2),vcup=vc(:,1),vcdn=vc(:,2))   