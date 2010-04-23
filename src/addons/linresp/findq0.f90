subroutine findq0(qvec,q0,a0)
use modmain
implicit none
real(8), intent(in) :: qvec(3,3)
real(8), intent(out) :: q0(3)
real(8), intent(out) :: a0
integer i1,i2,i3,n
real(8) vol,x0,x1,y0,y1,z0,z1,t
real(8), external :: r3mdet
real(8) v000(3),v001(3),v010(3),v100(3),v011(3),v101(3),v110(3),v111(3)
real(8) v1(3)

vol=abs(r3mdet(qvec))
n=100
t=0.d0
do i1=0,n-1
  x0=dble(i1)/n
  x1=dble(i1+1)/n
  do i2=0,n-1
    y0=dble(i2)/n
    y1=dble(i2+1)/n
    do i3=0,n-1
      z0=dble(i3)/n
      z1=dble(i3+1)/n
      v000(:)=x0*qvec(:,1)+y0*qvec(:,2)+z0*qvec(:,3)
      v001(:)=x0*qvec(:,1)+y0*qvec(:,2)+z1*qvec(:,3)
      v010(:)=x0*qvec(:,1)+y1*qvec(:,2)+z0*qvec(:,3)
      v100(:)=x1*qvec(:,1)+y0*qvec(:,2)+z0*qvec(:,3)
      v011(:)=x0*qvec(:,1)+y1*qvec(:,2)+z1*qvec(:,3)
      v101(:)=x1*qvec(:,1)+y0*qvec(:,2)+z1*qvec(:,3)
      v110(:)=x1*qvec(:,1)+y1*qvec(:,2)+z0*qvec(:,3)
      v111(:)=x1*qvec(:,1)+y1*qvec(:,2)+z1*qvec(:,3)
      if (.not.(i1.eq.0.and.i2.eq.0.and.i3.eq.0)) then      
        t=t+0.125d0*(1.d0/dot_product(v000,v000)+&
          1.d0/dot_product(v001,v001)+1.d0/dot_product(v010,v010)+&
          1.d0/dot_product(v100,v100)+1.d0/dot_product(v011,v011)+&
          1.d0/dot_product(v101,v101)+1.d0/dot_product(v110,v110)+&
          1.d0/dot_product(v111,v111))
      endif
    enddo
  enddo
enddo
write(*,*)'qvec=',qvec
write(*,*)'vol=',vol
t=t/(n**3)
write(*,*)'t=',t
q0(:)=0.5d0*(qvec(:,1)+qvec(:,2)+qvec(:,3))
! (2Pi)^-3 \int_dq f(q) = V*f(q0)*a0
! (2Pi)^-3 \int_dq f(q) = (2Pi)^-3 \sum f(q_i) delta_q
! delta_q = V/N
! (2Pi)^-3 \sum f(q_i) V/N = V*f(q0)*a0
! a0= (2Pi)^-3 \sum f(q_i) /N/f(q0)
a0=dot_product(q0,q0)*t/(twopi**3)

write(*,*)'f(q_i)*V*a0=',a0*vol/dot_product(q0,q0)
write(*,*)'q0=',q0,'a0=',a0

return
end
