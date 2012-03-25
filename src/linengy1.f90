
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
real(8), allocatable :: u(:)
!
allocate(u(nx))
call srrse_bound_state(sol,zn,n,l,nx,x,v,enu,u)
deallocate(u)
return
end subroutine

subroutine srrse_bound_state(sol,zn,n,l,nx,x,v,enu,p)
use mod_util
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
real(8), intent(out) :: p(nx)
!
real(8), allocatable :: ve(:),cve(:,:),mp(:),cmp(:,:),q(:),q1(:)
real(8), allocatable :: rhor2(:)
real(8) denu,t1
integer s,sp,i,j,j1,iter,nn
!
allocate(ve(nx))
allocate(cve(4,nx))
allocate(mp(nx),cmp(4,nx))
allocate(q(nx),q1(nx))
mp=0.d0
cmp=0.d0
allocate(rhor2(nx))
! electron part of the potential
do i=1,nx
  ve(i)=v(i)-zn/x(i)
enddo
! approximate ve with cubic spline
call spline(nx,x,1,ve,cve)
denu=0.01d0
!
do iter=1,1000
  call srrse_integrate(.false.,.true.,sol,zn,l,nx,x,v,ve,cve,mp,cmp,enu,p,q,q1,nn)
  sp=s
  if ((nn.gt.(n-l-1))) then
    s=-1
  else
    s=1
  endif
  denu=s*abs(denu)
  if (s.ne.sp) then
    denu=denu*0.5d0
  else
    denu=denu*1.25d0
  endif
  if (abs(denu).lt.1d-10) goto 10
  enu=enu+denu
enddo
write(*,'("Warning(srrse_bound_state): energy is not found")')
10 continue
! find the turning point
do i=1,nx
  if (v(i).gt.enu) exit
enddo
! find the minimum value of the function starting from the turning point
t1=1d100
do j=i,nx
  if (abs(p(j)).lt.t1) then
    t1=abs(p(j))
    j1=j
  endif
enddo
p(j1:)=0.d0
rhor2(:)=p(:)*p(:)
t1=rintegrate(nx,x,rhor2,m=0)
p(:)=p(:)/sqrt(t1)
nn=0
do i=1,nx-1
  if (p(i)*p(i+1).lt.0.d0) nn=nn+1
enddo
if (nn.ne.(n-l-1)) then
  write(*,'("Warning(srrse_bound_state): wrong number of nodes")')
  write(*,'("  n : ",I1,"  l : ",I1,"  nn : ",I1)')n,l,nn
endif
deallocate(ve,cve,mp,cmp,rhor2,q,q1)
end subroutine srrse_bound_state

subroutine srrse_solve(sol,zn,enu,l,m,nx,x,v,p,hp)
implicit none
real(8), intent(in) :: sol
real(8), intent(in) :: zn
real(8), intent(in) :: enu
integer, intent(in) :: l
integer, intent(in) :: m
integer, intent(in) :: nx
real(8), intent(in) :: x(nx)
real(8), intent(in) :: v(nx)
real(8), intent(out) :: p(nx)
real(8), intent(out) :: hp(nx)
!
integer i,j,nn
real(8) t1,alph2
real(8), allocatable :: ve(:),cve(:,:)
real(8), allocatable :: mp(:),cmp(:,:),q(:),q1(:)
!
allocate(ve(nx))
allocate(cve(4,nx))
! electron part of the potential
do i=1,nx
  ve(i)=v(i)-zn/x(i)
enddo
! approximate ve with cubic spline
call spline(nx,x,1,ve,cve)
allocate(mp(nx))
allocate(cmp(4,nx))
allocate(q(nx))
allocate(q1(nx))
do j=0,m
  if (j.eq.0) then
    mp(:)=0.d0
    cmp(:,:)=0.d0
  else
    mp(:)=dble(j)*p(:)
    call spline(nx,x,1,mp,cmp)
  endif
  call srrse_integrate(.false.,.false.,sol,zn,l,nx,x,v,ve,cve,mp,cmp,enu,p,q,q1,nn)
enddo
call fderiv(1,nx,x,q,mp,cmp)
if (.false.) then
  alph2=1.d0/(sol**2)
else
  alph2=0.d0
endif
do i=1,nx
  t1=2.d0-v(i)*alph2
  hp(i)=(dble(l*(l+1))/t1/(x(i)**2)+v(i))*p(i)-q(i)/x(i)-mp(i) !-q1(i)
enddo
deallocate(ve,cve,mp,cmp,q,q1)
return
end subroutine srrse_solve


!
! integrate scalar-relativistic radial Schrodinger equation using 4-th order 
!  Runge-Kutta method
!
subroutine srrse_integrate(trel,tfullm,sol,zn,l,nx,x,v,ve,cve,mp,cmp,enu,p,q,dqdr,nn)
implicit none
logical, intent(in) :: trel
logical, intent(in) :: tfullm
real(8), intent(in) :: sol
real(8), intent(in) :: zn
integer, intent(in) :: l
integer, intent(in) :: nx
real(8), intent(in) :: x(nx)
real(8), intent(in) :: v(nx)
real(8), intent(in) :: ve(nx)
real(8), intent(in) :: cve(4,nx)
real(8), intent(in) :: mp(nx)
real(8), intent(in) :: cmp(4,nx)
real(8), intent(in) :: enu
real(8), intent(out) :: p(nx)
real(8), intent(out) :: q(nx)
real(8), intent(out) :: dqdr(nx)
integer, intent(out) :: nn
!
integer i
real(8) k(2,4),h,h2,alph2,ll2,enu0
real(8) x0,x1,x2,m0,m1,m2,v0,v1,v2,t0,t1,t2,mp0,mp1,mp2,p0,p2,q0,q2
!
alph2=0.5d0*(1.d0/sol)**2
if (.not.trel) alph2=0.d0
enu0=0.d0
if (tfullm.and.trel) enu0=enu

ll2 = 0.5d0*dble(l*(l+1))

x2 = x(1)
v2 = v(1)
m2 = 1 - (v2-enu0)*alph2

p(1) = x(1)**(l+1)*exp(zn*x(1)/(l+1))
q(1) = (0.5d0/m2)*p(1)*(l/x(1)+zn/(l+1))

p2 = p(1)
q2 = q(1)
mp2 = mp(1)
t2 = ll2/m2/(x2**2)

do i=1,nx-1
  x0=x2
  x2=x(i+1)
  h=x2-x0
  h2=0.5d0*h
  x1=x0+h2
  p0=p2
  q0=q2
  m0=m2
  t0=t2
  v0=v2
  v2=v(i+1)
  mp0 = mp2
  mp2 = mp(i+1)
  mp1 = mp0 + cmp(1,i)*h2 + cmp(2,i)*(h2**2) + cmp(3,i)*(h2**3) 
  ! potential at the point x_i+h/2
  v1 = zn/x1 + ve(i) + cve(1,i)*h2 + cve(2,i)*(h2**2) + cve(3,i)*(h2**3)
  ! mass at the point x_i+h/2
  m1 = 1 - (v1-enu0)*alph2
  ! mass at the point x_i+h
  m2 = 1 - (v2-enu0)*alph2
  ! k1=F(Y(x),x)
  k(1,1) = 2*m0*q0 + p0/x0
  k(2,1) = (v0-enu+t0)*p0 - q0/x0 - mp0
  ! k2=F(Y(x)+k1*h/2,x+h/2)
  t1 = ll2/m1/(x1**2)
  k(1,2) = 2*m1*(q0+k(2,1)*h2) + (p0+k(1,1)*h2)/x1
  k(2,2) = (v1-enu+t1)*(p0+k(1,1)*h2) - (q0+k(2,1)*h2)/x1 - mp1
  ! k3=F(Y(x)+k2*h/2,x+h/2)
  k(1,3) = 2*m1*(q0+k(2,2)*h2) + (p0+k(1,2)*h2)/x1
  k(2,3) = (v1-enu+t1)*(p0+k(1,2)*h2) - (q0+k(2,2)*h2)/x1 - mp1
  ! k4=F(Y(x)+k3*h,x+h)
  t2=ll2/m2/(x2**2)
  k(1,4) = 2*m2*(q0+k(2,3)*h) + (p0+k(1,3)*h)/x2
  k(2,4) = (v2-enu+t2)*(p0+k(1,3)*h) - (q0+k(2,3)*h)/x2 - mp2
  ! Y(x+h)=Y(x)+h*(k1+2*k2+2*k3+k4)/6
  p2=p0+(k(1,1)+2*k(1,2)+2*k(1,3)+k(1,4))*h/6.d0
  dqdr(i)=(k(2,1)+2*k(2,2)+2*k(2,3)+k(2,4))/6.d0
  q2=q0+dqdr(i)*h
  p(i+1)=p2
  q(i+1)=q2
enddo
dqdr(nx)=dqdr(nx-1)
! count number of nodes
nn=0
do i=1,nx-1
  if (p(i)*p(i+1).lt.0.d0) nn=nn+1
enddo
return
end subroutine srrse_integrate


