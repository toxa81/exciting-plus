subroutine sic_genrmesh
use modmain
use mod_sic
implicit none
integer nrpole,npt
real(8), allocatable :: rpole(:)
integer maxdens
integer nrpt,i,j,ir,is
real(8) r0,x,alpha,d1,r0ratio,x0,a,b,c,dens0,dens1
real(8), allocatable :: rpt(:,:)
integer, parameter :: imesh=1
real(8) x1,x2,x3

if (allocated(s_r)) deallocate(s_r)

if (imesh.eq.1) then
  x0=1d-7
  dens0=500.d0
  dens1=50.d0
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
  allocate(s_r(s_nr))
  x0=1d-9
  b=log(sic_wan_cutoff/x0)/(s_nr-1)
  a=x0/exp(b)
  do ir=1,s_nr
    s_r(ir)=a*exp(ir*b)
  enddo
endif
if (imesh.eq.3) then
  allocate(s_r(s_nr))
  x0=1d-9
  do ir=1,s_nr
    s_r(ir)=x0+(sic_wan_cutoff-x0)*dble(ir-1)/dble(s_nr-1)
  enddo
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
 
!do ir=1,s_nr-1
!  s_rw(ir)=s_rw(ir)+0.5d0*(s_r(ir+1)-s_r(ir))*s_r(ir)**2
!  s_rw(ir+1)=s_rw(ir+1)+0.5d0*(s_r(ir+1)-s_r(ir))*s_r(ir+1)**2
!enddo

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

