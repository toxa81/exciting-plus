subroutine sic_genrmesh
use modmain
use mod_sic
implicit none
integer nrpole,npt
real(8), allocatable :: rpole(:)
integer maxdens
integer nrpt,i,j,ir
real(8) r0,x,alpha,d1,r0ratio,x0,a,b
real(8), allocatable :: rpt(:,:)
integer, parameter :: imesh=1

if (allocated(s_r)) deallocate(s_r)

if (imesh.eq.1) then
  nrpole=1
  allocate(rpole(nrpole))
  rpole(1)=0.d0
!  rpole(2)=5.9182d0
!  rpole(3)=8.3696d0
!  rpole(4)=10.25062d0
!  rpole(5)=11.83640d0
!  rpole(6)=13.23349d0
  maxdens=700
  nrpt=int(2*sic_wan_cutoff*maxdens)

  r0=1.8d0
  r0ratio=1/3.d0
  alpha=-log(r0ratio)/r0

  x=0.d0
  allocate(rpt(2,nrpt))
  do i=1,nrpt
    x=dble(i-1)*sic_wan_cutoff/dble(nrpt-1)
    rpt(1,i)=x
    d1=0.d0
    do j=1,nrpole
      d1=d1+maxdens*exp(-alpha*abs(x-rpole(j)))
    enddo
    rpt(2,i)=d1
  enddo
  d1=0.d0
  do i=1,nrpt-1
    d1=d1+0.5d0*(rpt(2,i)+rpt(2,i+1))*sic_wan_cutoff/dble(nrpt-1)
  enddo
  s_nr=int(d1)
  if (allocated(s_r)) deallocate(s_r)
  allocate(s_r(s_nr))
  d1=0.d0
  j=0
  do i=1,nrpt-1
    d1=d1+0.5d0*(rpt(2,i)+rpt(2,i+1))*sic_wan_cutoff/dble(nrpt-1)
    if (int(d1).gt.j) then
      j=j+1
      if (j.gt.s_nr) then
        write(*,'("Error(sic_genrmesh): j.gt.s_nr")')
        call pstop
      endif
      s_r(j)=rpt(1,i)
    endif
  enddo
  if (j.ne.s_nr) then
    write(*,'("Error(sic_genrmesh): j.ne.s_nr")')
    call pstop
  endif
  if (mpi_grid_root()) then
    write(*,*)
    write(*,'("[sic_genrmesh] number of radial points : ",I6)')s_nr
  endif
endif
if (imesh.eq.2) then
  allocate(s_r(s_nr))
  x0=1d-6
  b=log(sic_wan_cutoff/x0)/(s_nr-1)
  a=x0/exp(b)
  do ir=1,s_nr
    s_r(ir)=a*exp(ir*b)
  enddo
endif
if (imesh.eq.3) then
  allocate(s_r(s_nr))
  x0=1d-6
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
do ir=1,s_nr
  write(99,*)ir,s_r(ir)
enddo
do ir=1,s_nr
  write(100,*)s_r(ir),s_rw(ir)
enddo
return
end

