subroutine sic_gensmesh
use modmain
use mod_sic
implicit none
integer itp,lm
real(8) a,tp(2)
real(8), allocatable :: glw(:),glx(:)
integer nt,np,it,ip,l,m,lm1
real(8) t1
complex(8) zt1
complex(8), allocatable :: clm(:)
!
nt=lmaxwan+1
np=2*lmaxwan+1
s_ntp=nt*np
if (allocated(s_tp)) deallocate(s_tp)
allocate(s_tp(2,s_ntp))
if (allocated(s_x)) deallocate(s_x)
allocate(s_x(3,s_ntp))
if (allocated(s_tpw)) deallocate(s_tpw)
allocate(s_tpw(s_ntp))
if (allocated(s_rlmf)) deallocate(s_rlmf)
allocate(s_rlmf(lmmaxwan,s_ntp))
if (allocated(s_ylmf)) deallocate(s_ylmf)
allocate(s_ylmf(lmmaxwan,s_ntp))
if (allocated(s_rlmb)) deallocate(s_rlmb)
allocate(s_rlmb(s_ntp,lmmaxwan))
if (allocated(s_ylmb)) deallocate(s_ylmb)
allocate(s_ylmb(s_ntp,lmmaxwan))
! spherical grid with Gauss quadrature
allocate(glw(nt),glx(nt))
call gaulegf(-1.d0,1.d0,glx,glw,nt)
itp=0
do it=1,nt
  do ip=1,np
    itp=itp+1
    s_tp(1,itp)=acos(glx(it))
    s_tp(2,itp)=twopi*(ip-1)/np
    s_tpw(itp)=twopi*glw(it)/np
    s_x(:,itp)=(/sin(s_tp(1,itp))*cos(s_tp(2,itp)),&
                 sin(s_tp(1,itp))*sin(s_tp(2,itp)),&
                 glx(it)/)
  enddo
enddo
deallocate(glx,glw)
! generate spherical harmonics
do itp=1,s_ntp
  call genrlm(lmaxwan,s_tp(1,itp),s_rlmf(1,itp))
  call genylm(lmaxwan,s_tp(1,itp),s_ylmf(1,itp))
  do lm=1,lmmaxwan
    s_rlmb(itp,lm)=s_rlmf(lm,itp)*s_tpw(itp)
    s_ylmb(itp,lm)=dconjg(s_ylmf(lm,itp))*s_tpw(itp) 
  enddo
enddo
! generate Lebedev mesh for muffin-tins
! TODO: change to Gauss
if (allocated(mt_spx)) deallocate(mt_spx)
allocate(mt_spx(3,mt_ntp))
if (allocated(mt_tpw)) deallocate(mt_tpw)
allocate(mt_tpw(mt_ntp))
if (allocated(mt_ylmf)) deallocate(mt_ylmf)
allocate(mt_ylmf(lmmaxvr,mt_ntp))
call leblaik(mt_ntp,mt_spx,mt_tpw)
do itp=1,mt_ntp
  mt_tpw(itp)=mt_tpw(itp)*fourpi
  call sphcrd(mt_spx(:,itp),a,tp)
  call genylm(lmaxvr,tp,mt_ylmf(1,itp)) 
enddo

!call gensmesh_lebedev
!call gensmesh_uniform
call gensmesh_cubed
!call gensmesh_symcrys
!call gensmesh_healpix
!call gensmesh_icos

!do l=0,lmaxwan
!  t1=0.d0
!  do m=-l,l
!    lm=idxlm(l,m)
!    zt1=zzero
!    do itp=1,s_ntp
!      zt1=zt1+s_ylmb(itp,lm)
!    enddo
!    t1=t1+abs(zt1)**2
!  enddo
!  write(500,*)l,t1
!enddo

!call bstop

!allocate(clm(lmmaxwan))
!do lm=1,lmmaxwan
!  clm=zzero
!  do lm1=1,lmmaxwan
!    do itp=1,s_ntp
!      clm(lm1)=clm(lm1)+s_ylmb(itp,lm1)*s_ylmf(lm,itp)
!    enddo
!  enddo
!  a=0.d0
!  do lm1=1,lmmaxwan
!    if (lm1.eq.lm) clm(lm1)=clm(lm1)-zone
!    a=a+abs(clm(lm1))
!  enddo
!  write(*,*)lm,a
!enddo
!deallocate(clm)
!call bstop
!call sic_gensmesh_healpix
!call sic_gensmesh_ico
!call nfsft_init(lmaxwan,s_ntp,stp)
return
end








subroutine sic_gensmesh_ico
use modmain
use mod_sic
implicit none
integer nf,nx,ne,i,itp,lm,lm1,k,l,m
real(8) a,x_(3),p1(3),p2(3),p3(3),t1,t2
integer, allocatable :: facets(:,:)
integer, allocatable :: edges(:,:)
real(8), allocatable :: x(:,:)
real(8), allocatable :: x1(:,:)
logical tfound
complex(8), allocatable :: clm(:)
complex(8) zt1
a=(1.d0+sqrt(5.d0))/2.d0
nx=12
allocate(x(3,nx))
x(:,1)= (/ 0.d0, 1.d0, a/)
x(:,2)= (/ 0.d0, 1.d0,-a/)
x(:,3)= (/ 0.d0,-1.d0, a/)
x(:,4)= (/ 0.d0,-1.d0,-a/)
x(:,5)= (/ 1.d0, a, 0.d0/)
x(:,6)= (/ 1.d0,-a, 0.d0/)
x(:,7)= (/-1.d0, a, 0.d0/)
x(:,8)= (/-1.d0,-a, 0.d0/)
x(:,9)= (/ a, 0.d0, 1.d0/)
x(:,10)=(/ a, 0.d0,-1.d0/)
x(:,11)=(/-a, 0.d0, 1.d0/)
x(:,12)=(/-a, 0.d0,-1.d0/)
do i=1,nx
  x(:,i)=x(:,i)/sqrt(sum(x(:,i)**2))
enddo

! icosahedron recursive division: 
!  number of vertices 10*4^k + 2
!  number of facets   20*4^k
!  number of edges    30*4^k
do k=0,3 !if (allocated(xk)) deallocate(xk)
  !allocate(xk(3,10*4**k+2))
  !xk(:,:)=x(:,:)
  if (allocated(facets)) deallocate(facets)
  allocate(facets(3,20*4**k))
  call find_facets(nx,x,nf,facets)
  if (allocated(edges)) deallocate(edges)
  allocate(edges(2,30*4**k))
  call find_edges(nf,facets,ne,edges)

  do i=1,nf
    p1(:)=x(:,facets(3,i))-x(:,facets(1,i))
    p2(:)=x(:,facets(3,i))-x(:,facets(2,i))
    call r3cross(p1,p2,p3)
    write(*,*)"face : ",i,"surf=",sqrt(dot_product(p3,p3))
  enddo

  !do i=1,1
  !  p1=x(:,facets(1,i))
  !  p2=x(:,facets(2,i))
  !  p3=x(:,facets(3,i))
  !  call facet_weight(p1,p2,p3)
  !enddo
  !call bstop

  
  if (allocated(x1)) deallocate(x1)
  allocate(x1(3,10*4**(k+1)+2))
  x1(:,1:nx)=x(:,:)
  do i=1,ne
    x_(:)=0.5d0*(x(:,edges(1,i))+x(:,edges(2,i)))
    x1(:,nx+i)=x_(:)/sqrt(sum(x_(:)**2))
  enddo
  if (allocated(x)) deallocate(x)
  allocate(x(3,10*4**(k+1)+2))
  x=x1
  nx=10*4**(k+1)+2
enddo

write(*,*)"nx=",nx
s_ntp=nx
if (allocated(s_tp)) deallocate(s_tp)
allocate(s_tp(2,s_ntp))
if (allocated(s_x)) deallocate(s_x)
allocate(s_x(3,s_ntp))
if (allocated(s_tpw)) deallocate(s_tpw)
allocate(s_tpw(s_ntp))
if (allocated(s_rlmf)) deallocate(s_rlmf)
allocate(s_rlmf(lmmaxwan,s_ntp))
if (allocated(s_ylmf)) deallocate(s_ylmf)
allocate(s_ylmf(lmmaxwan,s_ntp))
if (allocated(s_rlmb)) deallocate(s_rlmb)
allocate(s_rlmb(s_ntp,lmmaxwan))
if (allocated(s_ylmb)) deallocate(s_ylmb)
allocate(s_ylmb(s_ntp,lmmaxwan))
s_x(:,1:nx)=x(:,1:nx)
do itp=1,s_ntp
  call sphcrd(s_x(1,itp),a,s_tp(1,itp))
  s_tpw(itp)=fourpi/s_ntp
enddo

do itp=1,s_ntp
  call genrlm(lmaxwan,s_tp(1,itp),s_rlmf(1,itp))
  call genylm(lmaxwan,s_tp(1,itp),s_ylmf(1,itp))
  do lm=1,lmmaxwan
    s_rlmb(itp,lm)=s_rlmf(lm,itp)*s_tpw(itp)
    s_ylmb(itp,lm)=dconjg(s_ylmf(lm,itp))*s_tpw(itp) 
  enddo
enddo

do l=0,lmaxwan
  t1=0.d0
  do m=-l,l
    lm=idxlm(l,m)
    zt1=zzero
    do itp=1,s_ntp
      zt1=zt1+s_ylmb(itp,lm)
    enddo
    t1=t1+abs(zt1)**2
  enddo
  write(500,*)l,t1
enddo
call bstop



!allocate(clm(lmmaxwan))
!do lm=1,lmmaxwan
!  clm=zzero
!  do lm1=1,lmmaxwan
!    do itp=1,s_ntp
!      clm(lm1)=clm(lm1)+s_ylmb(itp,lm1)*s_ylmf(lm,itp)
!    enddo
!  enddo
!  a=0.d0
!  do lm1=1,lmmaxwan
!    if (lm1.eq.lm) clm(lm1)=clm(lm1)-zone
!    a=a+abs(clm(lm1))
!  enddo
!  write(*,*)lm,a
!enddo
!deallocate(clm)
call bstop



return
end subroutine







subroutine facet_weight(x1,x2,x3)
use modmain
implicit none
real(8), intent(in) :: x1(3)
real(8), intent(in) :: x2(3)
real(8), intent(in) :: x3(3)

real(8) r1(3),r2(3),r3(3),rt(3),m(3,3),tp(2),theta,phi,a,rv(3),v1(3),v2(3)
real(8) p2,t3,p3,r
integer i,j,n
real(8) b,t1,w

real(8) sinp2,cosp2,sint3,cost3,sinp3,cosp3
real(8) jdet

!real(8), external :: jdet

r1=x1; r2=x2; r3=x3
! -phi_1 rotation around z 
call sphcrd(r1,a,tp)
phi=tp(2)
m=0.d0
m(1,1)=cos(phi); m(2,2)=m(1,1); m(3,3)=1.d0
m(2,1)=-sin(phi); m(1,2)=-m(2,1)
call r3mv(m,r1,rt); r1=rt
call r3mv(m,r2,rt); r2=rt
call r3mv(m,r3,rt); r3=rt
! -theta_1 rotation around y
call sphcrd(r1,a,tp)
theta=tp(1)
m=0.d0
m(1,1)=cos(theta); m(3,3)=m(1,1); m(2,2)=1.d0
m(3,1)=sin(theta); m(1,3)=-m(3,1)
call r3mv(m,r1,rt); r1=rt
call r3mv(m,r2,rt); r2=rt
call r3mv(m,r3,rt); r3=rt
! -phi_2 rotation around z
call sphcrd(r2,a,tp)
phi=tp(2)
m=0.d0
m(1,1)=cos(phi); m(2,2)=m(1,1); m(3,3)=1.d0
m(2,1)=-sin(phi); m(1,2)=-m(2,1)
call r3mv(m,r1,rt); r1=rt
call r3mv(m,r2,rt); r2=rt
call r3mv(m,r3,rt); r3=rt
! Pi/2 rotation around y, then Pi/2 rotation around x
m=0.d0
m(1,3)=1.d0
m(2,1)=1.d0
m(3,2)=1.d0
call r3mv(m,r1,rt); r1=rt
call r3mv(m,r2,rt); r2=rt
call r3mv(m,r3,rt); r3=rt

write(*,*)"r1=",r1
write(*,*)"r2=",r2
write(*,*)"r3=",r3

    v1(:)=r1-r2
    v2(:)=r1-r3
    call r3cross(v1,v2,rv)
    write(*,*)"surf=",sqrt(dot_product(rv,rv))
 

call sphcrd(r2,a,tp)
p2=tp(2)
call sphcrd(r3,a,tp)
t3=tp(1)
p3=tp(2)


sinp2=sin(p2)
cosp2=cos(p2)
sint3=sin(t3)
cost3=cos(t3)
sinp3=sin(p3)
cosp3=cos(p3)

t1=0.d0
n=40000
do i=1,n
  a=dble((i-1))/(n-1)
  do j=1,i
    w=1.d0
    if (i.eq.1.and.j.eq.1) w=0.125d0
    if (i.eq.n.and.j.eq.n) w=0.125d0
    if (i.eq.n.and.j.eq.1) w=0.25d0
    if (i.gt.1.and.i.lt.n.and.j.eq.1) w=0.5d0
    if (i.eq.n.and.j.gt.1.and.j.lt.n) w=0.5d0
    if (j.eq.i.and.i.gt.1.and.i.lt.n) w=0.5d0
    b=dble(j-1)/(n-1)
    rv(:)=r1(:)+(r2(:)-r1(:))*a+(r3(:)-r2(:))*b
    call sphcrd(rv,r,tp)
    theta=tp(1)
    phi=tp(2)
    
    
    jdet=(cost3*sinp2)/((b**2*cost3**2+(1.d0-a+(a-b)*cosp2+b*cosp3*sint3)**2+& 
      ((a-b)*sinp2+b*sinp3*sint3)**2)**1.5*sqrt(1.d0-(b**2*cost3**2)/ &
      (b**2*cost3**2+(1.d0-a+(a-b)*cosp2+b*cosp3*sint3)**2+((a-b)*sinp2+b*sinp3*sint3)**2)))
    
    t1=t1+sin(theta)*abs(jdet)*w/n/n
  enddo
enddo
write(*,*)"t1=",t1,fourpi/20




return
end subroutine


real(8) function jdet(a,b,p2,t3,p3)
implicit none
real(8), intent(in) :: a
real(8), intent(in) :: b
real(8), intent(in) :: p2
real(8), intent(in) :: t3
real(8), intent(in) :: p3
real(8) t
t=(Cos(t3)*Sin(p2))/ &
    ((b**2*Cos(t3)**2 + (1 - a + (a - b)*Cos(p2) + b*Cos(p3)*Sin(t3))**2 + & 
    ((a - b)*Sin(p2) + b*Sin(p3)*Sin(t3))**2)**1.5* &
     Sqrt(1 - (b**2*Cos(t3)**2)/ &
      (b**2*Cos(t3)**2 + (1 - a + (a - b)*Cos(p2) + b*Cos(p3)*Sin(t3))**2 +& 
           ((a - b)*Sin(p2) + b*Sin(p3)*Sin(t3))**2)))
jdet=t
return
end function


subroutine find_edges(nf,f,ne,e)
use modmain
implicit none
integer, intent(in) :: nf
integer, intent(in) :: f(3,nf)
integer, intent(out) :: ne
integer, intent(inout) :: e(2,*)
integer i,j,ie,ed(2,3),k
logical tfound

ie=0
do i=1,nf
  ed(1,1)=f(1,i)
  ed(2,1)=f(2,i)
  ed(1,2)=f(2,i)
  ed(2,2)=f(3,i)
  ed(1,3)=f(3,i)
  ed(2,3)=f(1,i)
  do k=1,3
    tfound=.false.
    do j=1,ie
      if ((e(1,j).eq.ed(1,k).and.e(2,j).eq.ed(2,k)).or.&
          (e(2,j).eq.ed(1,k).and.e(1,j).eq.ed(2,k)))  tfound=.true.
    enddo
    if (.not.tfound) then
      ie=ie+1
      e(:,ie)=ed(:,k)
    endif
  enddo
enddo
ne=ie
write(*,*)"number of edges : ",ie
return
end subroutine

subroutine find_facets(nx,x,nf,f)
use modmain
implicit none
integer, intent(in) :: nx
real(8), intent(in) :: x(3,nx)
integer, intent(out) :: nf
integer, intent(inout) :: f(3,*)
integer ifct,i1,i2,i3,n1,n2,j,j1,j2,nfp,facet(3)
real(8) x0(3),x1(3),x2(3),x3(3),norm(3),prod,mprod
logical tfound
integer, allocatable :: facet_(:)
allocate(facet_(nx))
ifct=0
do i1=1,nx-2
  do i2=i1+1,nx-1
    do i3=i2+1,nx
      x1(:)=x(:,i1)-x(:,i3)
      x2(:)=x(:,i2)-x(:,i3)
      call r3cross(x1,x2,norm)
      x0=0.d0
      j1=nx+1
      n1=0
      n2=0
      nfp=0
      do j=1,nx
        x3(:)=x(:,j)-x(:,i3)
        prod=dot_product(x3,norm)
        if (abs(prod).lt.1d-10) then
          nfp=nfp+1
          facet_(nfp)=j
          if (j.lt.j1) j1=j
          x0(:)=x0(:)+x(:,j)
        else if (prod.ge.1d-10) then
          n1=n1+1
        else 
          n2=n2+1
        endif
      enddo !j
      if (n1.eq.0.or.n2.eq.0) then
        if (nfp.ne.3) then
          stop "nfp.ne.3"
        endif
        if (n2.eq.0)  norm=-1.d0*norm
        x0=x0/nfp
        do n1=1,nfp
          x1(:)=x(:,j1)-x0(:)
          mprod=-100.d0
          do n2=1,nfp
            x2(:)=x(:,facet_(n2))-x0(:)
            call r3cross(x1,x2,x3)
            prod=dot_product(x1,x2)
            if (prod.gt.mprod.and.dot_product(x3,norm).gt.0.d0) then
              mprod=prod
              j2=facet_(n2)
            endif
          enddo !n2
          facet(n1)=j1
          j1=j2
        enddo !n1
        tfound=.false.
        do j=1,ifct
          if (all(f(:,j).eq.facet(:))) tfound=.true.
        enddo
        if (.not.tfound) then
          ifct=ifct+1
          f(:,ifct)=facet(:)
        endif
      endif
    enddo !i3
  enddo !i2
enddo !i1
write(*,*)"number of facets : ",ifct
nf=ifct
deallocate(facet_)
return
end subroutine


!subroutine sic_gensmesh_sym
!use modmain
!use mod_sic
!implicit none
!integer i,j,j1,j2,i1,i2,i3,nx,ias,jas,nf,n1,n2,nfp,itp
!integer facet_(100),facet(100)
!logical lfound
!real(8) x0(3),x1(3),x2(3),x3(3),norm(3),prod,mprod,a
!real(8) t1,p1,t2,p2,t3,p3,s
!real(8), allocatable :: x(:,:)
!integer, allocatable :: facets(:,:)
!logical tfound
!
!s_ntp=5
!
!call getnghbr(-0.d0,50.d0)
!allocate(x(3,nnghbr(1)))
!x=0.d0
!nx=0
!ias=1
!i1=inghbr(6,nnghbr(ias),ias)
!do j=2,nnghbr(ias)
!  if (inghbr(6,j,ias).lt.i1) then
!    jas=inghbr(1,j,ias)
!    x1(:)=atposc(:,ias2ia(jas),ias2is(jas))+&
!      inghbr(3,j,ias)*avec(:,1)+&
!      inghbr(4,j,ias)*avec(:,2)+&
!      inghbr(5,j,ias)*avec(:,3)-atposc(:,ias2ia(ias),ias2is(ias))
!    x1=x1/sqrt(sum(x1**2))
!    lfound=.false.
!    do i=1,nx
!      if (sum(abs(x1(:)-x(:,i))).lt.1d-10) lfound=.true.
!    enddo
!    if (.not.lfound) then
!      nx=nx+1
!      x(:,nx)=x1(:)
!    endif
!    if (nx.ge.s_ntp.and.inghbr(6,j+1,ias).gt.inghbr(6,j,ias)) exit
!  endif
!enddo
!if (nx.lt.s_ntp) then
!  write(*,'("Error(sic_gensmesh): not enough covering points")')
!  write(*,'("  found : ",I6)')nx
!  write(*,'("  minimum : ",I6)')s_ntp
!  call pstop
!endif
!if (mpi_grid_root()) then
!  write(*,'("[sic_gensmesh] number of covering points : ",I4)')nx
!endif
!
!allocate(facets(0:100,nx))
!nf=0
!do i1=1,nx-2
!  do i2=i1+1,nx-1
!    do i3=i2+1,nx
!      x1(:)=x(:,i1)-x(:,i3)
!      x2(:)=x(:,i2)-x(:,i3)
!      call r3cross(x1,x2,norm)
!      x0=0.d0
!      j1=nx+1
!      n1=0
!      n2=0
!      nfp=0
!      do j=1,nx
!        x3(:)=x(:,j)-x(:,i3)
!        prod=dot_product(x3,norm)
!        if (abs(prod).lt.1d-10) then
!          nfp=nfp+1
!          facet_(nfp)=j
!          if (j.lt.j1) j1=j
!          x0(:)=x0(:)+x(:,j)
!        else if (prod.ge.1d-10) then
!          n1=n1+1
!        else 
!          n2=n2+1
!        endif
!      enddo !j
!      if (n1.eq.0.or.n2.eq.0) then
!        if (n2.eq.0)  norm=-1.d0*norm
!        x0=x0/nfp
!        do n1=1,nfp
!          x1(:)=x(:,j1)-x0(:)
!          mprod=-100.d0
!          do n2=1,nfp
!            x2(:)=x(:,facet_(n2))-x0(:)
!            call r3cross(x1,x2,x3)
!            prod=dot_product(x1,x2)
!            if (prod.gt.mprod.and.dot_product(x3,norm).gt.0.d0) then
!              mprod=prod
!              j2=facet_(n2)
!            endif
!          enddo !n2
!          facet(n1)=j1
!          j1=j2
!        enddo !n1
!        tfound=.false.
!        do j=1,nf
!          if (facets(0,j).eq.nfp) then
!            if (all(facets(1:nfp,j).eq.facet(1:nfp))) tfound=.true.
!          endif
!        enddo
!        if (.not.tfound) then
!          nf=nf+1
!          facets(0,nf)=nfp
!          facets(1:nfp,nf)=facet(1:nfp)
!        endif
!      endif
!    enddo !i3
!  enddo !i2
!enddo !i1
!write(*,*)"number of facets : ",nf
!if (any(facets(0,1:nf).ne.3)) stop "facet is not triangle"
!s_ntp=nx
!
!
!
!
!
!
!
!
!if (allocated(s_tp)) deallocate(s_tp)
!allocate(s_tp(2,s_ntp))
!if (allocated(s_x)) deallocate(s_x)
!allocate(s_x(3,s_ntp))
!if (allocated(s_tpw)) deallocate(s_tpw)
!allocate(s_tpw(s_ntp))
!!if (allocated(s_rlmf)) deallocate(s_rlmf)
!!allocate(s_rlmf(lmmaxwan,s_ntp))
!!if (allocated(s_ylmf)) deallocate(s_ylmf)
!!allocate(s_ylmf(lmmaxwan,s_ntp))
!!if (allocated(s_rlmb)) deallocate(s_rlmb)
!!allocate(s_rlmb(s_ntp,lmmaxwan))
!!if (allocated(s_ylmb)) deallocate(s_ylmb)
!!allocate(s_ylmb(s_ntp,lmmaxwan))
!!
!do itp=1,s_ntp
!  s_x(:,itp)=x(:,itp)
!  call sphcrd(s_x(1,itp),a,s_tp(1,itp))
!enddo
!
!s_tpw=0.d0
!
!do j=1,nf
!  i1=facets(1,j)
!  i2=facets(2,j)
!  i3=facets(3,j)
!  t1=s_tp(1,i1)
!  p1=s_tp(2,i1)
!  t2=s_tp(1,i2)
!  p2=s_tp(2,i2)
!  t3=s_tp(1,i3)
!  p3=s_tp(2,i3)
!
!  !write(*,*)"before sort:",t1,t2,t3
!  if (t1.gt.t3) then
!    a=t1
!    t1=t3
!    t3=a
!    a=p1
!    p1=p3
!    p3=a
!    i=i1
!    i1=i3
!    i3=i
!  endif
!  if (t2.gt.t3) then
!    a=t2
!    t2=t3
!    t3=a
!    a=p2
!    p2=p3
!    p3=a
!    i=i2
!    i2=i3
!    i3=i
!  endif
!  if (t1.gt.t2) then
!    a=t1
!    t1=t2
!    t2=a
!    a=p1
!    p1=p2
!    p2=a
!    i=i1
!    i1=i2
!    i2=i
!  endif
!  write(*,*)"after sort:",t1,t2,t3,p1,p2,p3
!
!  s=abs(p2*t1 - p3*t1 - p1*t2 + p3*t2 + p1*t3 - p2*t3)
!  if (abs(t1-t2).lt.1d-10) then
!    s_tpw(i1)=s_tpw(i1)+s*(-((-2 + t2**2 - 2*t2*t3 + t3**2)*Cos(t2) + 2*(Cos(t3) + (-t2 + t3)*Sin(t2)))/(2.*(t2 - t3)**3))
!    s_tpw(i2)=s_tpw(i2)+s*(-((-2 + t2**2 - 2*t2*t3 + t3**2)*Cos(t2) + 2*(Cos(t3) + (-t2 + t3)*Sin(t2)))/(2.*(t2 - t3)**3))
!    s_tpw(i3)=s_tpw(i3)+s*((-2*Cos(t2) + 2*Cos(t3) - (t2 - t3)*(Sin(t2) + Sin(t3)))/(t2 - t3)**3)
!  else if (abs(t2-t3).lt.1d-10) then
!    s_tpw(i1)=s_tpw(i1)+s*((-2*Cos(t1) + 2*Cos(t3) - (t1 - t3)*(Sin(t1) + Sin(t3)))/(t1 - t3)**3)
!    s_tpw(i2)=s_tpw(i2)+s*((2*Cos(t1) + (-2 + t1**2 - 2*t1*t3 + t3**2)*Cos(t3) + 2*(t1 - t3)*Sin(t3))/(2.*(t1 - t3)**3))  
!    s_tpw(i3)=s_tpw(i3)+s*((2*Cos(t1) + (-2 + t1**2 - 2*t1*t3 + t3**2)*Cos(t3) + 2*(t1 - t3)*Sin(t3))/(2.*(t1 - t3)**3))
!  else
!    s_tpw(i1)=s_tpw(i1)+s*(((-2*t1*t2 + t2**2 + 2*t1*t3 - t3**2)*Cos(t1) + (t1 - t3)**2*Cos(t2) - &
!      (t1 - t2)*((t1 - t2)*Cos(t3) + (t1 - t3)*(t2 - t3)*Sin(t1)))/((t1 - t2)**2*(t1 - t3)**2*(t2 - t3)))
!    s_tpw(i2)=s_tpw(i2)+s*(((t2 - t3)**2*Cos(t1) + (t1 - t3)*(t1 - 2*t2 + t3)*Cos(t2) + (t1 - t2)*((-t1 + t2)*Cos(t3) + &
!      (t1 - t3)*(t2 - t3)*Sin(t2)))/((t1 - t2)**2*(t1 - t3)*(t2 - t3)**2))
!    s_tpw(i3)=s_tpw(i3)+s*(((t2 - t3)**2*Cos(t1) - (t1 - t3)**2*Cos(t2) + (t1 - t2)*((t1 + t2 - 2*t3)*Cos(t3) + &
!      (t1 - t3)*(-t2 + t3)*Sin(t3)))/((t1 - t2)*(t1 - t3)**2*(t2 - t3)**2))
!  endif
!enddo
!
!write(*,*)s_tpw
!write(*,*)sum(s_tpw)
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!
!call bstop
!return
!end
