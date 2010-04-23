subroutine init_q0gamma
use modmain
implicit none

integer i1,i2,i3,n
real(8) vq0l_1(3),vq0l_2(3),vq0l_3(3)
real(8) vq0c_1(3),vq0c_2(3),vq0c_3(3)
real(8) q0mtrx(3,3),vq0c__1(3),vq0c__2(3)

vq0l_1=(/1,0,0/)/dble(ngridk(1))
call r3mv(bvec,vq0l_1,vq0c_1)
vq0l_2=(/0,1,0/)/dble(ngridk(2))
call r3mv(bvec,vq0l_2,vq0c_2)
vq0l_3=(/0,0,1/)/dble(ngridk(3))
call r3mv(bvec,vq0l_3,vq0c_3)

n=0
do i1=0,1
  do i2=0,1
    do i3=0,1
      q0mtrx(:,1)=(i1-0.5d0)*vq0c_1(:)
      q0mtrx(:,2)=(i2-0.5d0)*vq0c_2(:)
      q0mtrx(:,3)=(i3-0.5d0)*vq0c_3(:)
      n=n+1
      call findq0(q0mtrx,q0gamma(1,n),a0gamma(n))
    enddo
  enddo
enddo
return
end