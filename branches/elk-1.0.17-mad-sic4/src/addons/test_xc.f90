subroutine test_xc
use modmain
implicit none
real(8) dens(nspinor),vxc(nspinor),exc
integer i
call init0

do i=0,10000
  dens=0.d0
  dens(1)=dble(i)*0.01
  call elk_xc(dens,vxc,exc)
  write(100,*)dens(1),exc
enddo

return
end
