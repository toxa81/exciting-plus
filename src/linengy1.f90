
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: linengy
! !INTERFACE:
subroutine linengy1
! !USES:
use modmain
! !DESCRIPTION:
!   Calculates the new linearisation energies for both the APW and local-orbital
!   radial functions. See the routine {\tt findband}.
!
! !REVISION HISTORY:
!   Created May 2003 (JKD)
!EOP
!BOC
implicit none
! local variables
logical fnd
integer is,ia,ja,ias,jas,ic,ir
integer l,ilo,io,jo
real(8) t1
real(8), allocatable :: vr(:)
e0min=100.d0
allocate(vr(spnrmax))
! loop over non-equivalent atoms (classes)
do ic=1,natmcls
  ias=ic2ias(ic)
  ia=ias2ia(ias)
  is=ias2is(ias)
! potential
  vr(1:nrmt(is))=veffmt(1,1:nrmt(is),ias)*y00
  t1=vr(nrmt(is))-spvr(nrmt(is),is)
  do ir=nrmt(is)+1,spnr(is)
    vr(ir)=spvr(ir,is)+t1
  enddo
  do l=0,lmaxapw
    do io=1,apword(l,is)
      if (apwve(io,l,is)) then
        call findband2(solsc,spzn(is),apwpqn(l,is),l,spnr(is),spr(1,is),vr,apwe(io,l,ias))
        e0min=min(e0min,apwe(io,l,ias))
      endif
    enddo
  enddo
  do ilo=1,nlorb(is)
    do io=1,lorbord(ilo,is)
      if (lorbve(io,ilo,is)) then
        l=lorbl(ilo,is)
        call findband2(solsc,spzn(is),lopqn(ilo,is),l,spnr(is),spr(1,is),vr,lorbe(io,ilo,ias))
        e0min=min(e0min,lorbe(io,ilo,ias))
      endif
    enddo
  enddo
  do jas=1,natmtot
    if (ias2ic(jas).eq.ic) then
      do l=0,lmaxapw
        do io=1,apword(l,is)
          apwe(io,l,jas)=apwe(io,l,ias)
        enddo
      enddo
      do ilo=1,nlorb(is)
        do io=1,lorbord(ilo,is)
          lorbe(io,ilo,jas)=lorbe(io,ilo,ias)
        enddo
      enddo
    endif
  enddo
enddo
deallocate(vr)
return
end subroutine
!EOC



subroutine findband2(sol,zn,n,l,nx,x,v,enu)
implicit none
real(8), intent(in) :: sol
real(8), intent(in) :: zn
integer, intent(in) :: n
integer, intent(in) :: l
integer, intent(in) :: nx
real(8), intent(in) :: x(nx)
real(8), intent(in) :: v(nx)
real(8), intent(out) :: enu
!
real(8), allocatable :: u(:),up(:)
!
allocate(u(nx),up(nx))
call bound_state(sol,zn,n,l,nx,x,v,enu,u,up)
deallocate(u,up)
return
end subroutine

subroutine bound_state(sol,zn,n,l,nx,x,v,enu,y,yp)
use mod_sic
implicit none
!
real(8), intent(in) :: sol
real(8), intent(in) :: zn
integer, intent(in) :: n
integer, intent(in) :: l
integer, intent(in) :: nx
real(8), intent(in) :: x(nx)
real(8), intent(in) :: v(nx)
real(8), intent(out) :: enu
real(8), intent(out) :: y(nx)
real(8), intent(out) :: yp(nx)
!
real(8), allocatable :: v1(:),cv1(:,:)
real(8), allocatable :: rhor2(:)
real(8) denu,t1
integer s,sp,i,j,j1,iter,nn
!
allocate(v1(nx))
allocate(cv1(4,nx))
allocate(rhor2(nx))
! unknown part of the potential
do i=1,nx
  v1(i)=v(i)-zn/x(i)
enddo
! approximate v1 with cubic spline
call spline(nx,x,1,v1,cv1)
enu=v(nx)-1.d0
denu=0.5d0
!
do iter=1,2000
  call integrate_srrseq_rk4(sol,zn,l,nx,x,v,v1,cv1,enu,y,yp,nn)
  sp=s
  if ((nn.gt.(n-l-1))) then
    s=-1
  else
    s=1
  endif
  denu=s*abs(denu)
  if (s.ne.sp) then
    denu=denu*0.5d0
  endif
  if (abs(denu).lt.1d-16) goto 10
  enu=enu+denu
enddo
write(*,'("Warning(bound_state): energy is not found")')
10 continue
! find the turning point
do i=1,nx
  if (v(i).gt.enu) exit
enddo
! find the minimum value of the function starting from the turning point
t1=1d100
do j=i,nx
  if (abs(y(j)).lt.t1) then
    t1=abs(y(j))
    j1=j
  endif
enddo
y(j1:)=0.d0
yp(j1:)=0.d0
rhor2(:)=y(:)*y(:)
t1=rintegrate(nx,x,rhor2,m=0)
y(:)=y(:)/sqrt(t1)
yp(:)=yp(:)/sqrt(t1)
nn=0
do i=1,nx-1
  if (y(i)*y(i+1).lt.0.d0) nn=nn+1
enddo
deallocate(v1,cv1,rhor2)
end subroutine bound_state

!
! integrate scalar-relativistic radial Schrodinger equation using 4-th order 
!  Runge-Kutta method
!
subroutine integrate_srrseq_rk4(sol,zn,l,nx,x,v,v1,cv1,enu,u,up,nn)
implicit none
real(8), intent(in) :: sol
real(8), intent(in) :: zn
integer, intent(in) :: l
integer, intent(in) :: nx
real(8), intent(in) :: x(nx)
real(8), intent(in) :: v(nx)
real(8), intent(in) :: v1(nx)
real(8), intent(in) :: cv1(4,nx)
real(8), intent(in) :: enu
real(8), intent(out) :: u(nx)
real(8), intent(out) :: up(nx)
integer, intent(out) :: nn
integer i
real(8), allocatable :: rm(:),drm(:),vl(:)
real(8) h,h2,m1,m2,m3,m4,k1,k2,k3,k4,v2,v2p
real(8) q0,q1,q2,p0,p1,p2,a,b,alph2

u(1)=x(1)**(l+1)*exp(zn*x(1)/(l+1))
up(1)=(l+1)*x(1)**l*exp(zn*x(1)/(l+1))+zn*u(1)/(l+1)

alph2=0.5d0*(1.d0/sol)**2
alph2=alph2/1000.d0

allocate(rm(nx),drm(nx),vl(nx))
do i=1,nx
  rm(i)=1.d0+(enu-v(i))*alph2
  drm(i)=-alph2*(-zn/(x(i)**2)+cv1(1,i))
  vl(i)=dble(l*(l+1))/(x(i)**2)
enddo

do i=1,nx-1
  h=x(i+1)-x(i)
  h2=h/2.d0
  ! "mass" function at the point x_i
  q0=rm(i)
  ! potential at the point x_i + h/2
  v2=zn/(x(i)+h2)+v1(i)+cv1(1,i)*h2+cv1(2,i)*(h2**2)+cv1(3,i)*(h2**3)
  ! "mass" function at the point x_i + h/2
  q1=1.d0+(enu-v2)*alph2
  ! "mass" function at the point x_{i+1}
  q2=rm(i+1)
  ! "mass" function derivative at the point x_i
  p0=drm(i)
  ! potential derivative at the point x_i + h/2
  v2p=-zn/((x(i)+h2)**2)+cv1(1,i)+2.d0*cv1(2,i)*h2+3.d0*cv1(3,i)*(h2**2)
  ! "mass" function derivative at the point x_i + h/2
  p1=-alph2*v2p
  ! "mass" function derivative at the point x_{i+1}
  p2=drm(i+1)
  !
  m1=up(i)
  b=p0/q0
  a=2.d0*q0*(v(i)-enu)+vl(i)-b/x(i)
  k1=a*u(i)+b*up(i)
  !
  m2=up(i)+h2*k1
  b=p1/q1
  a=2.d0*q1*(v2-enu)+dble(l*(l+1))/((x(i)+h2)**2)-b/(x(i)+h2)
  k2=a*(u(i)+h2*m1)+b*(up(i)+h2*k1)
  !
  m3=up(i)+h2*k2
  k3=a*(u(i)+h2*m2)+b*(up(i)+h2*k2)
  !
  m4=up(i)+h*k3
  b=p2/q2
  a=2.d0*q2*(v(i+1)-enu)+vl(i+1)-b/x(i+1)
  k4=a*(u(i)+h*m3)+b*(up(i)+h*k3)
  !
  u(i+1)=u(i)+(m1+2*m2+2*m3+m4)*h/6.d0
  up(i+1)=up(i)+(k1+2*k2+2*k3+k4)*h/6.d0
enddo
nn=0
do i=1,nx-1
  if (u(i)*u(i+1).lt.0.d0) nn=nn+1
enddo
deallocate(rm,drm,vl)
return
end subroutine


