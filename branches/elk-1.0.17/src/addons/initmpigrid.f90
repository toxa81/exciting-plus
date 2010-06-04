subroutine initmpigrid
use modmain
use mod_addons_q
implicit none
integer, allocatable :: sz(:)
integer i1,i2
!------------------------!
!     parallel grid      !
!------------------------!
if (task.eq.0.or.task.eq.1.or.task.eq.22) then
  i2=2
  allocate(sz(i2))
  if (nproc.le.nkpt) then
    sz=(/nproc,1/)
  else
    i1=nproc/nkpt
    sz=(/nkpt,i1/)
  endif    
else if (task.eq.400.or.task.eq.401.or.task.eq.802.or.task.eq.810) then
  i2=3
  allocate(sz(i2))
  if (nproc.le.nkptnr) then
    sz=(/nproc,1,1/)
  else
    i1=nproc/nkptnr
    if (i1.le.nvq) then
      sz=(/nkptnr,1,i1/)
    else
      sz=(/nkptnr,nproc/(nkptnr*nvq),nvq/)
    endif
  endif
else
  i2=1
  allocate(sz(i2))
  sz=(/nproc/)
endif  
! overwrite default grid layout
if (lmpigrid) then
  sz(1:i2)=mpigrid(1:i2) 
endif
call mpi_grid_initialize(sz)
deallocate(sz)
return
end