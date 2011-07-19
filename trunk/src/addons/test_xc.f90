subroutine test_xc
use modmain
implicit none
real(8) dens(nspinor),vx(nspinor),vc(nspinor),ex,ec
integer i
call init0

do i=0,10000
  dens=0.d0
  dens(1)=dble(i)*0.01
  call elk_xc(dens,vx,vc,ex,ec)
  write(100,*)dens(1),vx(1),ex
  write(200,*)dens(1),vc(1),ec
enddo

return
end
