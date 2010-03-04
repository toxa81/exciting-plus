subroutine init_band_trans
use modmain
implicit none
integer ikloc,ik,i,j,n1,n2,i1,i2,i3
real(8) t1
logical l1
logical, external :: wann_diel

! get limits
call getmeidx(.true.)
if (allocated(nmegqblh)) deallocate(nmegqblh)
allocate(nmegqblh(nkptnrloc))
nmegqblh=0
if (allocated(bmegqblh)) deallocate(bmegqblh)
allocate(bmegqblh(2,nmegqblhmax,nkptnrloc))
bmegqblh=0
if (wannier_megq) then
  if (allocated(nmegqblhwan)) deallocate(nmegqblhwan)
  allocate(nmegqblhwan(nkptnrloc))
  nmegqblhwan=0
  if (allocated(imegqblhwan)) deallocate(imegqblhwan)
  allocate(imegqblhwan(nmegqblhmax,nkptnrloc))
  imegqblhwan=0
endif
! setup n,n' stuff
call getmeidx(.false.)
! split interband transitions between second dimension
if (allocated(nmegqblhloc)) deallocate(nmegqblhloc)
allocate(nmegqblhloc(2,nkptnrloc))
nmegqblhloc=0
do ikloc=1,nkptnrloc
  nmegqblhloc(1,ikloc)=mpi_grid_map(nmegqblh(ikloc),dim_b,offs=i)
  nmegqblhloc(2,ikloc)=i
enddo
nmegqblhlocmax=maxval(nmegqblhloc)
if (wannier_megq) then
  nmegqblhwanmax=maxval(nmegqblhwan)
endif
return
end