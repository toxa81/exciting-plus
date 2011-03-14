subroutine sic_genrmesh
use modmain
use mod_sic
implicit none
integer npole,npt
real(8), allocatable :: rpole(:)
real(8), allocatable :: dpole(:)
integer maxdens
integer nrpt,i,j,ir,is,nr2,nr4
real(8) r0,x,alpha,d1,r0ratio,x0,a,b,c,dens0,dens1,t,dx,d,norm
real(8), allocatable :: rpt(:,:)
integer, parameter :: imesh=4
real(8) x1,x2,x3
real(8), external :: x_aux

if (allocated(s_r)) deallocate(s_r)

if (imesh.eq.1) then
  x0=1d-7
  dens0=5000.d0
  dens1=10.d0
  a=(dens1-dens0)/(sic_wan_cutoff-x0)
  b=dens1-a*sic_wan_cutoff
  c=1-x0*(0.5d0*a*x0+b)
  s_nr=int(sic_wan_cutoff*(0.5d0*a*sic_wan_cutoff+b)+c)+1
  allocate(s_r(s_nr))
  do ir=1,s_nr
    s_r(ir)=(-b+sqrt(b**2-2*a*c+2*a*ir))/a
  enddo
  s_r(1)=x0
  s_r(s_nr)=sic_wan_cutoff
endif
if (imesh.eq.2) then
  x0=1d-7
  nr2=s_nr/2
  nr4=nr2/2
! position of the nearest neigbour
  x1=8.36960d0 !5.91820d0
  x2=(x1-x0)/2.d0
  b=log(x2)-log(x0)
  allocate(s_r(s_nr))
  do ir=1,nr2
   t=(ir-1)/dble(nr2-1)
   s_r(ir)=x0*exp(b*t)
  enddo
  do ir=nr2+1,nr2+nr4
   t=(nr2+nr4-ir)/dble(nr4)
   s_r(ir)=x1-x0*exp(b*abs(t))
  enddo
  do ir=nr2+nr4+1,s_nr
    t=(ir-nr2-nr4)/dble(nr4)
    s_r(ir)=x1+x0*exp(b*abs(t))
  enddo
  sic_wan_cutoff=s_r(s_nr)
!  b=log(sic_wan_cutoff)-log(x0)
endif
if (imesh.eq.3) then
  allocate(s_r(s_nr))
  x0=1d-9
  do ir=1,s_nr
    s_r(ir)=x0+(sic_wan_cutoff-x0)*dble(ir-1)/dble(s_nr-1)
  enddo
endif
if (imesh.eq.4) then
  npole=2
  allocate(rpole(npole))
  rpole(1)=0.d0
  rpole(2)=5.91820d0
  !rpole(3)=8.36960d0
  !rpole(4)=10.25062d0
  allocate(s_r(s_nr))
  b=15.d0
  do ir=1,s_nr
    t=(dble(npole)*(ir-1)/dble(s_nr-1))-1.d0
    x=x_aux(1.d0*(npole-1),2.d0*(sic_wan_cutoff-rpole(npole)),b,t)
    do j=2,npole
      x=x+x_aux(1.d0*(j-2),rpole(j)-rpole(j-1),b,t)
    enddo
    s_r(ir)=x
  enddo
  deallocate(rpole)
  if (mpi_grid_root()) then
    write(*,'("[sic_genrmesh] first radial point : ",G18.10)')s_r(1)
  endif
endif
! generate radial weights for integration
if (allocated(s_rw)) deallocate(s_rw)
allocate(s_rw(s_nr)); s_rw=0.d0
x2=s_r(1)
x3=s_r(2)
s_rw(1)=-((x2-x3)*(3*x2**2+2*x2*x3+x3**2))/12.d0
do ir=2,s_nr-1
  x1=s_r(ir-1)
  x2=s_r(ir)
  x3=s_r(ir+1)
  s_rw(ir)=-((x1-x3)*(x1**2+x2**2+x2*x3+x3**2+x1*(x2+x3)))/12.d0
enddo
x1=s_r(s_nr-1)
x2=s_r(s_nr)
s_rw(s_nr)=-((x1-x2)*(x1**2+2*x1*x2+3*x2**2))/12.d0
! radial weights for muffin-tins
if (allocated(mt_rw)) deallocate(mt_rw)
allocate(mt_rw(nrmtmax,nspecies))
do is=1,nspecies
  x2=spr(1,is)
  x3=spr(2,is)
  mt_rw(1,is)=-((x2-x3)*(3*x2**2+2*x2*x3+x3**2))/12.d0
  do ir=2,nrmt(is)-1
    x1=spr(ir-1,is)
    x2=spr(ir,is)
    x3=spr(ir+1,is)
    mt_rw(ir,is)=-((x1-x3)*(x1**2+x2**2+x2*x3+x3**2+x1*(x2+x3)))/12.d0
  enddo
  x1=spr(nrmt(is)-1,is)
  x2=spr(nrmt(is),is)
  mt_rw(nrmt(is),is)=-((x1-x2)*(x1**2+2*x1*x2+3*x2**2))/12.d0
  if (sum(mt_rw(1:nrmt(is),is))-((fourpi/3)*rmt(is)**3).gt.1d-10) then
    write(*,'("Error(sic_genrmesh): wrong weight for is : ",I4)')is
    call pstop
  endif
enddo 
return
end

real(8) function x_aux(t0,h,b,t)
implicit none
real(8), intent(in) :: t0
real(8), intent(in) :: h
real(8), intent(in) :: b
real(8), intent(in) :: t
x_aux=h*(1.d0-1.d0/(exp(b*(t-t0))+1.d0))
return
end function

