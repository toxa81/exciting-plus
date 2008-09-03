subroutine wann_plot
use modmain
use modwann
implicit none

real(8) r(3),r1(3),t(3),d

integer ntr(3),i,ivec

call init0
call init1


r(:)=(/10.1d0,0.8973d0,0.445555d0/)
write(*,*)avec
write(*,*)bvec

do ivec=1,3
  d=dot_product(r,bvec(:,ivec))/twopi
  ntr(ivec)=floor(d)
enddo

t(:)=ntr(1)*avec(:,1)+ntr(2)*avec(:,2)+ntr(3)*avec(:,3)

r1(:)=r(:)-t(:)

write(*,*)'ntr=',ntr
write(*,*)'r=',r
write(*,*)'r1=',r1
write(*,*)'t=',t


return
end
