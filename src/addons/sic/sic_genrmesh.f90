subroutine sic_genrmesh
use modmain
use mod_sic
implicit none
integer is,i,j,ir
real(8) x,x0,x1,x2,x3,b,t
integer, parameter :: imesh=2
real(8), external :: x_aux

if (allocated(s_r)) deallocate(s_r)

! linear mesh
if (imesh.eq.1) then
  allocate(s_r(s_nr))
  x0=1d-9
  do ir=1,s_nr
    s_r(ir)=x0+(sic_wan_cutoff-x0)*dble(ir-1)/dble(s_nr-1)
  enddo
endif
! multi-pole mesh
if (imesh.eq.2) then
  if (.not.allocated(s_rpole)) then
    s_nrpole=1
    allocate(s_rpole(1))
    s_rpole=0.d0
  endif
  do i=1,s_nrpole
    if (s_rpole(i).ge.sic_wan_cutoff) then
      write(*,'("Error(sic_genrmesh): pole of the radial mesh is greater than cutoff radius")')
      write(*,'("  pole(",I1,") : ",G18.10)')i-1,s_rpole(i)
      write(*,'("  sic_wan_cutoff : ",G18.10)')sic_wan_cutoff
      call pstop
    endif
  enddo
  allocate(s_r(s_nr))
  b=15.d0
  do ir=1,s_nr
! t is in interval [-1, s_nrpole-1]
    t=(dble(s_nrpole)*(ir-1)/dble(s_nr-1))-1.d0
    x=x_aux(1.d0*(s_nrpole-1),2.d0*(sic_wan_cutoff-s_rpole(s_nrpole)),b,t)
    do j=2,s_nrpole
      x=x+x_aux(1.d0*(j-2),s_rpole(j)-s_rpole(j-1),b,t)
    enddo
    s_r(ir)=x
  enddo
  if (abs(s_r(s_nr)-sic_wan_cutoff).gt.1d-5) then
    write(*,'("Error(sic_genrmesh): last point of radial mesh is wrong")')
    write(*,'("  s_r(s_nr) : ",G18.10)')s_r(s_nr)
    write(*,'("  sic_wan_cutoff : ",G18.10)')sic_wan_cutoff
    call pstop
  endif
  s_r(s_nr)=sic_wan_cutoff
  if (mpi_grid_root()) then
    write(*,'("[sic_genrmesh] first radial point : ",G18.10)')s_r(1)
  endif
endif
! generate radial weights for integration
if (allocated(s_rw)) deallocate(s_rw)
allocate(s_rw(s_nr))
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

