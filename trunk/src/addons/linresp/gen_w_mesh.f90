subroutine gen_w_mesh
use modmain
use mod_linresp
use mod_addons
implicit none
integer iw
!
if (lr_nw.eq.1) then
  lr_dw=0.d0
else
  if (timgw) then
    lr_dw=(lr_iw1-lr_iw0)/(lr_nw-1)
  else
    lr_dw=(lr_w1-lr_w0)/(lr_nw-1)
  endif
endif
if (allocated(lr_w)) deallocate(lr_w)
allocate(lr_w(lr_nw))
do iw=1,lr_nw
  if (timgw) then
    lr_w(iw)=zi*(lr_iw0+lr_dw*(iw-1))
  else
    lr_w(iw)=dcmplx(lr_w0+lr_dw*(iw-1),lr_eta)
  endif
enddo
return
end subroutine
