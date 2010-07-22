subroutine init_q0
use modmain
use mod_addons_q
implicit none
integer i1,i2,i3,n
real(8) qvec(3,3)
n=0
do i1=0,1
  do i2=0,1
    do i3=0,1
      qvec(:,1)=(i1-0.5d0)*bvec(:,1)/ngridk(1)
      qvec(:,2)=(i2-0.5d0)*bvec(:,2)/ngridk(2)
      qvec(:,3)=(i3-0.5d0)*bvec(:,3)/ngridk(3)
      n=n+1
      call findq0(qvec,vq0c(1,n),aq0(n))
    enddo
  enddo
enddo
return
end