program qdir
use modmain
implicit none

real(8) a(3,3),b(3),q(3),q1(3),ang
integer i,j,k,lwork,ipiv(3)
integer nq(3)
real(8) work(24)

! q-vector, direction.
q=(/1.130773d0, 0.000000d0, -31.573183d0/)
! normalize q-vector
q(:)=q(:)/sqrt(sum(q(:)**2))
! length in inverse angstroms converted to inverse a.u.
q(:)=3.4*q(:)*au2ang

write(*,*)'q=',q

call readinput
call init0

do j=1,3
  do i=1,j
    a(i,j)=dot_product(bvec(:,i)/ngridk(i),bvec(:,j)/ngridk(j))
  enddo
  b(j)=dot_product(bvec(:,j)/ngridk(j),q(:))
enddo
lwork=24
call dsysv('U',3,1,a,3,ipiv,b,3,work,lwork,i)
write(*,*)'f1,f2,f3=',b
do i=1,3
  nq(i)=floor(b(i))
enddo

do i=-ngridk(1),ngridk(1)
do j=-ngridk(2),ngridk(2)
do k=-ngridk(3),ngridk(3)
  q1(:)=(nq(1)+i)*bvec(:,1)/ngridk(1)+&
        (nq(2)+j)*bvec(:,2)/ngridk(2)+&
        (nq(3)+k)*bvec(:,3)/ngridk(3)
  ang=acos((q1(1)*q(1)+q1(2)*q(2)+q1(3)*q(3))/sqrt(sum(q1**2))/sqrt(sum(q**2)))
  write(*,'("angle : ",F8.4," diff : ",F8.4," coord : ",3I4)')ang,sqrt(sum((q-q1)**2)),nq(1)+i,nq(2)+j,nq(3)+k
enddo
enddo
enddo
        
  




return
end
