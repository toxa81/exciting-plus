subroutine sic_test_mesh
use modmain
use mod_sic
implicit none
real(8), external :: d_trial_func
complex(8), external :: z_trial_func
real(8), external :: d_trial_radial_func
real(8), allocatable :: ftp(:,:),flm(:,:)
complex(8), allocatable :: zftp(:,:),zflm(:,:)
real(8) x(3),t1
integer itp,ir,lm

allocate(ftp(s_ntp,s_nr))
allocate(flm(lmmaxwan,s_nr))
allocate(zftp(s_ntp,s_nr))
allocate(zflm(lmmaxwan,s_nr))
do ir=1,s_nr
  do itp=1,s_ntp
    x(:)=s_spx(:,itp)*s_r(ir)
    zftp(itp,ir)=z_trial_func(x,0,1)
  enddo
enddo

!call dgemm('T','N',lmmaxwan,s_nr,s_ntp,1.d0,s_rlmb,s_ntp,ftp,&
!  s_ntp,0.d0,flm,lmmaxwan)

call zgemm('T','N',lmmaxwan,s_nr,s_ntp,zone,s_ylmb,s_ntp,zftp,&
  s_ntp,zzero,zflm,lmmaxwan)

t1=0.d0
do lm=1,lmmaxwan
  do ir=1,s_nr
    t1=t1+abs(zflm(lm,ir))**2*s_rw(ir)
  enddo
enddo
write(*,*)t1
write(*,*)s_dot_ll((/0.d0,0.d0,0.d0/),(/2.d0,0.d0,0.d0/),zflm,zflm) 

!t1=0.d0
!do ir=1,s_nr
!  t1=t1+abs(flm(17,ir)-d_trial_radial_func(s_r(ir),4))**2
!enddo
!write(*,*)"diff=",t1
!
!t1=0.d0
!do ir=1,s_nr
!  t1=t1+abs(flm(1,ir))**2*s_rw(ir)
!enddo
!write(*,*)t1


call pstop
return
end

real(8) function d_trial_func(x,l,m)
implicit none
real(8), intent(in) :: x(3)
integer, intent(in) :: l
integer, intent(in) :: m
real(8), external :: d_trial_radial_func
real(8) r,tp(2)
real(8) rlm(100)
call sphcrd(x,r,tp)
call genrlm(l,tp,rlm)
d_trial_func=d_trial_radial_func(r,l)*rlm(l**2+m)
return
end

real(8) function d_trial_radial_func(r,l)
use modmain
implicit none
real(8), intent(in) :: r
integer, intent(in) :: l
!d_trial_radial_func=(r**l)*exp(-2.d0*r)
d_trial_radial_func=exp(-0.5*r**2)*cos(r*twopi/5.d0)
return
end

complex(8) function z_trial_func(x,l,m)
use modmain
implicit none
real(8), intent(in) :: x(3)
integer, intent(in) :: l
integer, intent(in) :: m
real(8), external :: d_trial_radial_func
real(8) r,tp(2)
real(8) rlm(100)
call sphcrd(x,r,tp)
call genrlm(l,tp,rlm)
z_trial_func=zone*d_trial_radial_func(r,l) !*rlm(l**2+m)
return
end

