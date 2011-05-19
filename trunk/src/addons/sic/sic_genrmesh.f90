subroutine sic_genrmesh
use modmain
use mod_sic
implicit none
integer is,i,j,ir
real(8) x,x0,x1,x2,x3,t,a
real(8), allocatable :: b(:)
integer, parameter :: imesh=2
real(8), external :: x_aux

if (allocated(s_r)) deallocate(s_r)

! linear mesh
if (imesh.eq.1) then
  allocate(s_r(s_nr))
  x0=1d-9
  do ir=1,s_nr
    s_r(ir)=x0+(s_rmax-x0)*dble(ir-1)/dble(s_nr-1)
  enddo
endif
! multi-pole mesh
if (imesh.eq.2) then
  if (.not.allocated(s_rpole)) then
    s_nrpole=1
    allocate(s_rpole(1))
    s_rpole=0.d0
  endif
! find poles 
  call getnghbr(-0.d0,s_rmax-0.1d0)
  s_nrpole=inghbr(6,nnghbr(1),1)
  deallocate(s_rpole)
  allocate(s_rpole(s_nrpole))
  do i=1,s_nrpole
    do j=1,nnghbr(1)
      if (inghbr(6,j,1).eq.i) then
        s_rpole(i)=inghbr(2,j,1)/1000000.d0
        exit
      endif
    enddo !j
  enddo !i
  do i=1,s_nrpole
    if (s_rpole(i).ge.s_rmax) then
      write(*,'("Error(sic_genrmesh): pole of the radial mesh is greater &
        &than sphere radius")')
      write(*,'("  pole(",I1,") : ",G18.10)')i-1,s_rpole(i)
      write(*,'("  s_rmax : ",G18.10)')s_rmax
      call pstop
    endif
  enddo
  if (mpi_grid_root()) then
    write(*,'("[sic_genrmesh] poles of the radial mesh : ")')
    write(*,'(20F12.6)')s_rpole
  endif
  allocate(b(s_nrpole))
  b=15.d0
  !b(1)=16.d0
! compute a=x(s_nr)
  t=dble(s_nrpole)-1.d0
  a=x_aux(1.d0*(s_nrpole-1),2.d0*(s_rmax-s_rpole(s_nrpole)),&
    b(s_nrpole),t)
  do j=2,s_nrpole
    a=a+x_aux(1.d0*(j-2),s_rpole(j)-s_rpole(j-1),b(j-1),t)
  enddo
  allocate(s_r(s_nr))
  do ir=1,s_nr
! t is in interval [-1, s_nrpole-1]
    t=(dble(s_nrpole)*(ir-1)/dble(s_nr-1))-1.d0
    x=x_aux(1.d0*(s_nrpole-1),2.d0*(s_rmax-s_rpole(s_nrpole)),&
      b(s_nrpole),t)
    do j=2,s_nrpole
      x=x+x_aux(1.d0*(j-2),s_rpole(j)-s_rpole(j-1),b(j-1),t)
    enddo
    s_r(ir)=x*s_rmax/a
  enddo
  deallocate(b)
endif
! find number of radial points for s_rmin 
s_nr_min=s_nr
do ir=s_nr,1,-1
  if (s_r(ir).gt.s_rmin) then
    s_nr_min=ir
  endif
enddo
if (mpi_grid_root()) then
  write(*,'("[sic_genrmesh] first and last radial points : ")')
  write(*,'(2G18.10)')s_r(1),s_r(s_nr)
  write(*,'("[sic_genrmesh] number of points for s_rmin : ",I8)')s_nr_min
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

