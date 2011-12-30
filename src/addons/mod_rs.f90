!> @brief Radial solver
!! @details Generate radial solutions and radial integrals of the LAPW.\n
!! Preliminary notes:\n
!! 1. atomic wave-function \f$\psi({\bf r})\f$ is separated into radial and angular parts:
!!  \f$\psi({\bf r}) = R_{nl}(r)Y_{lm}(\theta,\phi)\f$ \n
!! 2. Laplacian in spherical coordinates: 
!!  \f$\Delta = \frac{1}{r^2} \frac{\partial}{\partial r} \big( r^2 \frac{\partial}{\partial r} \big) + \frac{1}{r^2} \Delta_{\theta,\phi}\f$\n
!! 3. Spherical harmonics are the eigenfunctions of the surface Laplacian: \f$ \Delta_{\theta,\phi}Y_{lm}(\theta,\phi) = -l(l+1)Y_{lm}(\theta,\phi)\f$\n
!! 4. Eigen-value equation: \f$ -\frac{1}{2}\big( \Delta_{r} + \frac{1}{r^2} \Delta_{\theta,\phi} \big)R_{nl}(r)Y_{lm}(\theta,\phi)+
!!  V(r)R_{nl}(r)Y_{lm}(\theta,\phi)=\epsilon R_{nl}(r)Y_{lm}(\theta,\phi) \f$\n
!! 5. Radial part of the equation: \f$ -\frac{1}{2} \frac{1}{r^2} \frac{\partial}{\partial r} \big( r^2 \frac{\partial R_{nl}(r)}{\partial r} \big)+
!!  \frac{l(l+1)}{2r^2}R_{nl}(r)+V(r)R_{nl}(r)=\epsilon R_{nl}(r)\f$
!! .
!! @author Anton Kozhevnikov
!! @date 2011
module mod_rs

!
! Notes 
! 1. atomic wave-function \psi(r) is separated into radial and 
!    angular parts: \psi(r) = R_{nl}(r)*Y_{lm}(\theta,\phi)
! 2. the following change is done: R=u/r which leads to the well-known
!    radial equation: -1/2 u"(r) + (V(r)+l*(l+1)/2/r^2)*u(r)=E*u(r)

!> number of mesh points
integer rs_nx
!> radial mesh
real(8), allocatable :: rs_x(:)
!> r^2 integration weights
real(8), allocatable :: rs_w(:)
!> trapezoidal integration weights
real(8), allocatable :: rs_wt(:)
contains

!> @brief Initialize radial grid and integration weights
subroutine rs_init
implicit none
real(8) x0,x1,x2,x3,t1
integer depth,k,j,n,i,n2

x0=1d-6
depth=16
n2=500

rs_nx=(n2*2+1)+depth*n2

allocate(rs_x(rs_nx))
allocate(rs_w(rs_nx))
allocate(rs_wt(rs_nx))
rs_x(1)=x0
k=1
do j=depth,0,-1
  if (j.eq.depth) then
    x1=x0
    x2=x0+(1.d0-x0)/(2**j)
    n=2*n2
  else
    x1=x0+(1.d0-x0)/(2**(j+1))
    x2=x0+(1.d0-x0)/(2**j)
    n=n2
  endif
  do i=1,n
    k=k+1
    rs_x(k)=x1+(x2-x1)*dble(i)/n
  enddo
enddo
do i=1,rs_nx-1
  if (rs_x(i).gt.rs_x(i+1)) then
    write(*,'("Error[rs_init]: mesh points are in the wrong order")')
    stop
  endif
enddo
x2=rs_x(1)
x3=rs_x(2)
rs_w(1)=-((x2-x3)*(3*x2**2+2*x2*x3+x3**2))/12.d0
rs_wt(1)=0.5d0*(x3-x2)
do i=2,rs_nx-1
  x1=rs_x(i-1)
  x2=rs_x(i)
  x3=rs_x(i+1)
  rs_w(i)=-((x1-x3)*(x1**2+x2**2+x2*x3+x3**2+x1*(x2+x3)))/12.d0
  rs_wt(i)=0.5d0*(x3-x1)
enddo
x1=rs_x(rs_nx-1)
x2=rs_x(rs_nx)
rs_w(rs_nx)=-((x1-x2)*(x1**2+2*x1*x2+3*x2**2))/12.d0
rs_wt(rs_nx)=0.5d0*(x2-x1)
t1=0.d0
do i=1,rs_nx
  t1=t1+rs_w(i)
enddo
if (abs(t1-(1.d0-x0**3)/3).gt.1d-14) then
  write(*,'("Error[rs_init]: wrong radial weights")')
  stop
endif
return
end subroutine

!> @brief Integrate second-order differential equation
subroutine rse_integrate(nr,r,f,g,u)
implicit none
integer, intent(in) :: nr
real(8), intent(in) :: r(nr)
real(8), intent(in) :: f(nr)
real(8), intent(in) :: g(nr)
real(8), intent(inout) :: u(nr)
!
integer i,i1,i2
real(8) h,h2,d1,d2
!
! solve u"(x)=f(x)u(x)+g(x)
! u(1) and u(2) must be initialized
do i=3,nr
  h=r(i)-r(i-1)
  i1=i-1
  i2=i-2
  if (abs(r(i1)-r(i2)-h).gt.1d-10) i2=i-3
  h2=h**2
  d1=-u(i2)+2*u(i1)+h2*(f(i1)*u(i1)+g(i1))+h2*(g(i)+f(i2)*u(i2)+g(i2)-2*f(i1)*u(i1)-2*g(i1))/12.d0
  d2=1-h2*f(i)/12.d0
  u(i)=d1/d2
enddo
return
end subroutine

!> @brief Integrate m-th energy derivative of the radial equation  
subroutine rs_solve(nr,l,m,zn,enu,r,v,u,nn)
implicit none
integer, intent(in) :: nr
integer, intent(in) :: l
integer, intent(in) :: m
real(8), intent(in) :: zn
real(8), intent(in) :: enu
real(8), intent(in) :: r(nr)
real(8), intent(in) :: v(nr)
real(8), intent(inout) :: u(nr)
integer, intent(out) :: nn
!
integer ir,j
real(8) t1,a,b
real(8), allocatable :: f(:),g(:)
!
allocate(f(nr))
allocate(g(nr))
do ir=1,nr
  f(ir)=2.d0*(v(ir)+0.5d0*l*(l+1)/r(ir)**2-enu)
enddo

do j=0,m
  if (j.eq.0) then
    g=0.d0
  else
    do ir=1,nr
      g(ir)=-2.d0*dble(j)*u(ir)
    enddo
  endif
  do ir=1,2
    a=-abs(zn)/(l+1)
    b=(-a**2-2*enu)/(4*l+6)
    u(ir)=exp(a*r(ir)+b*r(ir)**2)*r(ir)**(l+1)
  enddo
  call rse_integrate(nr,r,f,g,u)
  !t1=0.d0
  !do ir=1,nr
  !  t1=t1+(r(nr)**3)*rse_rw(ir)*(u(ir)/r(ir))**2
  !enddo
  !write(*,*)"norm=",t1
  !u(:)=u(:)/sqrt(t1)
enddo
! normalise
t1=0.d0
do ir=1,nr
  t1=t1+r(nr)*rs_wt(ir)*u(ir)**2
enddo
u(:)=u(:)/sqrt(t1)
! get number of nodes
nn=0
do ir=1,nr-1
  if (u(ir)*u(ir+1).lt.0.d0) nn=nn+1
enddo
deallocate(f,g)
end subroutine

!> @brief Compute <u1|H|u2>
subroutine rs_h0(l,x,vx,f1,f2,h0)
implicit none
integer, intent(in) :: l
real(8), intent(in) :: x(rs_nx)
real(8), intent(in) :: vx(rs_nx)
real(8), intent(in) :: f1(rs_nx)
real(8), intent(in) :: f2(rs_nx)
real(8), intent(out) :: h0
!
real(8), allocatable :: df1(:),df2(:),cf(:,:)
integer i
!
h0=0.d0
allocate(df1(rs_nx))
allocate(df2(rs_nx))
allocate(cf(4,rs_nx))

call fderiv(1,rs_nx,x,f1,df1,cf)
call fderiv(1,rs_nx,x,f2,df2,cf)
do i=1,rs_nx
  h0=h0+df1(i)*df2(i)*rs_wt(i)*x(rs_nx)
enddo
h0=-0.5d0*(f1(rs_nx)*df2(rs_nx)-f1(1)*df2(1)-h0)
do i=1,rs_nx
  h0=h0+(vx(i)+0.5d0*l*(l+1)/x(i)**2)*f1(i)*f2(i)*rs_wt(i)*x(rs_nx)
enddo
deallocate(df1,df2,cf)
return
end subroutine

subroutine rs_linen(pqn,l,zn,enu,x,vx)
implicit none
integer, intent(in) :: pqn
integer, intent(in) :: l
real(8), intent(in) :: zn
real(8), intent(inout) :: enu
real(8), intent(in) :: x(rs_nx)
real(8), intent(in) :: vx(rs_nx)
!
real(8) t1
real(8), allocatable :: u(:)
integer maxstp
logical ttop,tbot
integer ie,nn,nnd,nndp,i
real(8) de,t,tp,etop,ebot

allocate(u(rs_nx))

!do i=1,30
!  enu=-32.d0 + dble(i)/10
!  call rs_solve(rs_nx,l,0,zn,enu,x,vx,u,nn)
!  call rs_h0(l,x,vx,u,u,t1)
!  write(*,*)"enu=",enu," nn=",nn," h0=",t1,"deriv=",(u(rs_nx)/x(rs_nx)-u(rs_nx-1)/x(rs_nx-1))/(x(rs_nx)-x(rs_nx-1))
!  do j=1,rs_nx
!    write(400+i,*)x(j),u(j)/x(j)
!  enddo
!enddo
! local variables
! maximum number of steps
de=0.5d0
etop=3.d0
ttop=.false.
tbot=.false.
maxstp=2000
! find the top of the band
do i=1,maxstp
  call rs_solve(rs_nx,l,0,zn,etop,x,vx,u,nn) 
  nnd=nn-(pqn-l-1)
  if (nnd.gt.0) then
    etop=etop-de
  else
   etop=etop+de
  end if
  if (i.gt.1) then
    if ((nnd.ne.0).or.(nndp.ne.0)) then
      if (nnd*nndp.le.0) then
        de=de*0.5d0
      else
        de=de*1.1d0
      end if
    end if
  end if
  nndp=nnd
  if (de.lt.1d-9) then
    ttop=.true.
    exit
  endif
enddo
de=-0.1d0
ebot=etop+3*abs(de)
do i=1,maxstp
  ebot=ebot+de
  call rs_solve(rs_nx,l,0,zn,ebot,x,vx,u,nn)
  t=(u(rs_nx)/x(rs_nx)-u(rs_nx-1)/x(rs_nx-1))/(x(rs_nx)-x(rs_nx-1))
  if (i.gt.1) then
    if (t*tp.le.0.d0) then
      if (abs(de).lt.1d-9) then
        tbot=.true.
        exit
      endif
      de=-0.5*de
    endif
  endif
  tp=t
enddo
if (ttop.and.tbot) then
! set the band energy to the mid-point
  enu=(etop+ebot)/2.d0
  return
endif
!write(*,*)"enu=",enu
!call rs_solve(rs_nx,l,0,zn,enu,x,vx,u,nn)
!do i=1,rs_nx
!  write(400,*)x(i),u(i)/x(i)
!enddo
!write(*,*)"etop=",etop
!call rs_solve(rs_nx,l,0,zn,etop,x,vx,u,nn)
!do i=1,rs_nx
!  write(401,*)x(i),u(i)/x(i)
!enddo
!write(*,*)"ebot=",ebot
!call rs_solve(rs_nx,l,0,zn,ebot,x,vx,u,nn)
!do i=1,rs_nx
!  write(402,*)x(i),u(i)/x(i)
!enddo
!stop "rs_linen"
deallocate(u)
return
end subroutine

subroutine genufr
use modmain
implicit none
! local variables
integer is,ia,ja,ias,jas
integer nr,ir,nn,l,io,jo,ic,j,ilo,info,np
integer l1,io1,l2,io2,l3,m3,lm3,ilo1,ilo2
real(8) t1
! automatic arrays
logical done(natmmax)
real(8), allocatable :: x(:),vx(:),u(:,:,:),ulo(:,:,:)
real(8), allocatable :: flo(:,:)
integer, allocatable :: ipiv(:)
real(8), allocatable :: xa(:),ya(:)
real(8), allocatable :: a(:,:),b(:),c(:)
real(8) vr(nrmtmax)
real(8) r2(nrmtmax),fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)
real(8), external :: polynom
!
haa=0.d0
hloa=0.d0
hlolo=0.d0
oalo=0.d0
ololo=0.d0
np=max(maxlorbord+1,4)
allocate(x(rs_nx))
allocate(vx(rs_nx))
allocate(u(rs_nx,apwordmax,0:lmaxapw))
allocate(ulo(rs_nx,maxlorbord,maxlorb))
allocate(flo(rs_nx,maxlorb))
allocate(ipiv(np))
allocate(xa(np),ya(np),c(np))
allocate(a(np,np),b(np))
do ic=1,natmcls
  ias=ic2ias(ic)
  is=ias2is(ias)
  nr=nrmt(is)
! use spherical part of the potential
  vr(1:nr)=veffmt(1,1:nr,ias)*y00
! fine mesh 
  x(:)=rs_x(:)*rmt(is)
! map V(r) from coarse to fine mesh
  call rfinterp(nr,spr(1,is),1,vr(1),rs_nx,x,1,vx(1))
! generate u(r) for APW
  do l=0,lmaxapw
    do io=1,apword(l,is)
      if (apwve(io,l,is)) then
        call rs_linen(apwpqn(l,is),l,spzn(is),apwe(io,l,ias),x,vx)
      endif
! integrate the radial Schrodinger equation
      call  rs_solve(rs_nx,l,apwdm(io,l,is),spzn(is),apwe(io,l,ias),&
        x,vx,u(1,io,l),nn)
! subtract linear combination of previous vectors
      do jo=1,io-1
        t1=0.d0
        do ir=1,rs_nx
          t1=t1+u(ir,io,l)*u(ir,jo,l)*rs_wt(ir)*rmt(is)
        enddo
        u(:,io,l)=u(:,io,l)-t1*u(:,jo,l)
      enddo
! normalise radial functions
      t1=0.d0
      do ir=1,rs_nx
        t1=t1+rs_wt(ir)*rmt(is)*u(ir,io,l)**2
      enddo
      if (t1.lt.1.d-10) then
        write(*,*)
        write(*,'("Error(genufr): degenerate APW radial functions")')
        write(*,'(" for species ",I4)') is
        write(*,'(" atom ",I4)') ia
        write(*,'(" angular momentum ",I4)') l
        write(*,'(" and order ",I4)') io
        write(*,*)
        stop
      endif
      u(:,io,l)=u(:,io,l)/sqrt(t1)
! derivative at the muffin-tin surface: R=u/r -> R'=u'/r-u/r^2
      t1=x(rs_nx)-x(rs_nx-1)
      apwdfr(io,l,ias)=(u(rs_nx,io,l)-u(rs_nx-1,io,l))/t1/rmt(is)-&
        u(rs_nx,io,l)/(rmt(is)**2)
! map to original mesh and divide by r
      call rfinterp(rs_nx,x,1,u(1,io,l),nr,spr(1,is),1,apwfr(1,1,io,l,ias))
      do ir=1,nr
        apwfr(ir,1,io,l,ias)=apwfr(ir,1,io,l,ias)/spr(ir,is)
      enddo
   enddo !io
  enddo !l
! local orbitals
  do ilo=1,nlorb(is)
    l=lorbl(ilo,is)
    do jo=1,lorbord(ilo,is)
      if (lorbve(jo,ilo,is)) then
        call rs_linen(lopqn(ilo,is),l,spzn(is),lorbe(jo,ilo,ias),x,vx)
      endif
! integrate the radial Schrodinger equation
      call  rs_solve(rs_nx,l,lorbdm(jo,ilo,is),spzn(is),lorbe(jo,ilo,ias),&
        x,vx,ulo(1,jo,ilo),nn)
! set up the matrix of radial derivatives
      do j=1,np
        ir=rs_nx-np+j
        xa(j)=x(ir)
        ya(j)=ulo(ir,jo,ilo)/x(ir)
      end do
      do io=1,lorbord(ilo,is)
        a(io,jo)=polynom(io-1,np,xa,ya,c,rmt(is))
      end do
    end do
! set up the target vector
    b(:)=0.d0
    b(lorbord(ilo,is))=1.d0
    call dgesv(lorbord(ilo,is),1,a,np,ipiv,b,np,info)
    if (info.ne.0) then
      write(*,*)
      write(*,'("Error(genufr): degenerate local-orbital radial &
       &functions")')
      write(*,'(" for species ",I4)') is
      write(*,'(" atom ",I4)') ia
      write(*,'(" and local-orbital ",I4)') ilo
      write(*,'(" ZGESV returned INFO = ",I8)') info
      write(*,*)
      stop
    end if
! generate linear superposition of radial functions
    flo(:,ilo)=0.d0
    do io=1,lorbord(ilo,is)
      t1=b(io)
      flo(:,ilo)=flo(:,ilo)+t1*ulo(:,io,ilo)
    enddo
! normalise radial functions
    t1=0.d0
    do ir=1,rs_nx
      t1=t1+rs_wt(ir)*rmt(is)*flo(ir,ilo)**2
    enddo
    flo(:,ilo)=flo(:,ilo)/sqrt(t1)
! map to original mesh and divide by r
    call rfinterp(rs_nx,x,1,flo(1,ilo),nr,spr(1,is),1,lofr(1,1,ilo,ias))
    do ir=1,nr
      lofr(ir,1,ilo,ias)=lofr(ir,1,ilo,ias)/spr(ir,is)
    enddo
  enddo !ilo
! spherical part of APW-APW integrals
  do l=0,lmaxapw
    do io1=1,apword(l,is)
      do io2=1,apword(l,is)
        call rs_h0(l,x,vx,u(:,io1,l),u(:,io2,l),t1)
        haa(1,io1,l,io2,l,ias)=t1/y00
        !write(*,*)"my haa",haa(1,io1,l,io2,l,ias)
      enddo
    enddo
  enddo
! spherical part of lo-APW integrals 
  do ilo=1,nlorb(is)
    l=lorbl(ilo,is)
    do io2=1,apword(l,is)
      call rs_h0(l,x,vx,flo(:,ilo),u(:,io2,l),t1)
      hloa(1,ilo,io2,l,ias)=t1/y00
      !write(*,*)"my hloa",hloa(1,ilo,io2,l,ias)
      t1=0.d0
      do ir=1,rs_nx
        t1=t1+u(ir,io2,l)*flo(ir,ilo)*rs_wt(ir)*rmt(is)
      enddo
      oalo(io2,ilo,ias)=t1
    enddo
  enddo
! spherical part of lo-lo integrals
  do ilo1=1,nlorb(is)
    l1=lorbl(ilo1,is)
    do ilo2=1,nlorb(is)
      l2=lorbl(ilo2,is)
      if (l1.eq.l2) then
        call rs_h0(l1,x,vx,flo(:,ilo1),flo(:,ilo2),t1)
        hlolo(1,ilo1,ilo2,ias)=t1/y00
        !write(*,*)"my hlolo",hlolo(1,ilo1,ilo2,ias)
        t1=0.d0
        do ir=1,rs_nx
          t1=t1+flo(ir,ilo1)*flo(ir,ilo2)*rs_wt(ir)*rmt(is)
        enddo
        ololo(ilo1,ilo2,ias)=t1
     endif
    enddo
  enddo
! copy spherical part to the remaining atoms
  do jas=1,natmtot
    if (ias2ic(jas).eq.ic) then
      haa(1,:,:,:,:,jas)=haa(1,:,:,:,:,ias)
      hloa(1,:,:,:,jas)=hloa(1,:,:,:,ias)
      hlolo(1,:,:,jas)=hlolo(1,:,:,ias)
      oalo(:,:,jas)=oalo(:,:,ias)
      ololo(:,:,jas)=ololo(:,:,ias)
      apwfr(:,1,:,:,jas)=apwfr(:,1,:,:,ias)
      lofr(:,1,:,jas)=lofr(:,1,:,ias)
    endif
  enddo
end do !ic
deallocate(x)
deallocate(vx)
deallocate(u)
deallocate(ulo)
deallocate(flo)
deallocate(ipiv)
deallocate(xa,ya,c)
deallocate(a,b)
! non-spherical terms
do is=1,nspecies
  nr=nrmt(is)
  do ir=1,nr
    r2(ir)=spr(ir,is)**2
  end do
  do ia=1,natoms(is)
    ias=idxas(ia,is)
!---------------------------!
!     APW-APW integrals     !
!---------------------------!
    do l1=0,lmaxapw
      do io1=1,apword(l1,is)
        do l2=0,lmaxapw
          do io2=1,apword(l2,is)
            do l3=1,lmaxvr
              do m3=-l3,l3
                lm3=idxlm(l3,m3)
                do ir=1,nr
                  t1=apwfr(ir,1,io1,l1,ias)*apwfr(ir,1,io2,l2,ias)*r2(ir)
                  fr(ir)=t1*veffmt(lm3,ir,ias)
                end do
                call fderiv(-1,nr,spr(:,is),fr,gr,cf)
                haa(lm3,io1,l1,io2,l2,ias)=gr(nr)
              end do
            end do
          end do
        end do
      end do
    end do
!-------------------------------------!
!     local-orbital-APW integrals     !
!-------------------------------------!
    do ilo=1,nlorb(is)
      do l2=0,lmaxapw
        do io2=1,apword(l2,is)
          do l3=1,lmaxvr
            do m3=-l3,l3
              lm3=idxlm(l3,m3)
              do ir=1,nr
                t1=lofr(ir,1,ilo,ias)*apwfr(ir,1,io2,l2,ias)*r2(ir)
                fr(ir)=t1*veffmt(lm3,ir,ias)
              end do
              call fderiv(-1,nr,spr(:,is),fr,gr,cf)
              hloa(lm3,ilo,io2,l2,ias)=gr(nr)
            end do
          end do
        end do
      end do
    end do
!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
    do ilo1=1,nlorb(is)
      do ilo2=1,nlorb(is)
        l2=lorbl(ilo2,is)
        do l3=1,lmaxvr
          do m3=-l3,l3
            lm3=idxlm(l3,m3)
            do ir=1,nr
              t1=lofr(ir,1,ilo1,ias)*lofr(ir,1,ilo2,ias)*r2(ir)
              fr(ir)=t1*veffmt(lm3,ir,ias)
            end do
            call fderiv(-1,nr,spr(:,is),fr,gr,cf)
            hlolo(lm3,ilo1,ilo2,ias)=gr(nr)
          end do
        end do
      end do
    end do
! end loops over atoms and species
  end do
end do

!do l=0,lmaxapw
!  do ir=1,rs_nx
!    write(100,*)x(ir),u(ir,1,l)
!  enddo
!  write(100,*)
!enddo
!do ilo=1,nlorb(1)
!  do ir=1,rs_nx
!    write(200,*)x(ir),flo(ir,ilo)/x(ir)
!  enddo
!  write(200,*)
!enddo

!do l=0,lmaxapw
!  do ir=1,nrmt(1)
!    write(200,*)spr(ir,1),apwfr(ir,1,1,l,1)
!  enddo
!  write(200,*)
!enddo
!stop "genufr"

!do ilo=1,nlorb(1)
!  do ir=1,nrmt(1)
!    write(210,*)spr(ir,1),lofr(ir,1,ilo,1)
!  enddo
!  write(210,*)
!enddo

return
end subroutine

end module
