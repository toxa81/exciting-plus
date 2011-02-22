subroutine sic_genrmesh
use modmain
use mod_sic
implicit none
integer nrpole,npt
real(8), allocatable :: rpole(:)
integer maxdens
integer nrpt,i,j,ir
real(8) r0,x,alpha,d1,r0ratio,x0,a,b,c,dens0,dens1
real(8), allocatable :: rpt(:,:)
integer, parameter :: imesh=1

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
allocate(s_rw(s_nr))
s_rw=0.d0
do ir=1,s_nr-1
  s_rw(ir)=s_rw(ir)+0.5d0*(s_r(ir+1)-s_r(ir))*s_r(ir)**2
  s_rw(ir+1)=s_rw(ir+1)+0.5d0*(s_r(ir+1)-s_r(ir))*s_r(ir+1)**2
enddo
return
end

