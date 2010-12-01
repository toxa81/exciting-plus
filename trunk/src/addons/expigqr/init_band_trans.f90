subroutine init_band_trans
use modmain
use mod_linresp
implicit none
integer ikloc,i
logical, external :: wann_diel

! get limits
call getmeidx(.true.)
if (allocated(nmegqblhtot)) deallocate(nmegqblhtot)
allocate(nmegqblhtot(nkptnrloc))
nmegqblhtot=0
if (allocated(bmegqblh)) deallocate(bmegqblh)
allocate(bmegqblh(2,nmegqblhtotmax,nkptnrloc))
bmegqblh=0
if (wannier_megq) then
  if (allocated(nmegqblhwan)) deallocate(nmegqblhwan)
  allocate(nmegqblhwan(nkptnrloc))
  nmegqblhwan=0
  if (allocated(imegqblhwan)) deallocate(imegqblhwan)
  allocate(imegqblhwan(nmegqblhtotmax,nkptnrloc))
  imegqblhwan=0
endif
! setup n,n' stuff
call getmeidx(.false.)
! split interband transitions between second dimension
if (allocated(nmegqblhloc)) deallocate(nmegqblhloc)
allocate(nmegqblhloc(2,nkptnrloc))
nmegqblhloc=0
do ikloc=1,nkptnrloc
  nmegqblhloc(1,ikloc)=mpi_grid_map(nmegqblhtot(ikloc),dim_b,offs=i)
  nmegqblhloc(2,ikloc)=i
enddo
nmegqblhlocmax=maxval(nmegqblhloc(1,:))
if (wannier_megq) then
  nmegqblhwanmax=maxval(nmegqblhwan)
endif
return
end