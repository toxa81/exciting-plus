subroutine initmpigrid
use modmain
use mod_addons_q
implicit none
integer, allocatable :: d(:)
integer i1,i2,nd
!------------------------!
!     parallel grid      !
!------------------------!
if (task.eq.0.or.task.eq.1.or.task.eq.22..or.task.eq.20.or.&
    task.eq.805.or.task.eq.807.or.task.eq.808.or.task.eq.811) then
  nd=2
  allocate(d(nd))
  d=1
  if (nproc.le.nkpt) then
    d(dim1)=nproc
  else
    d(dim1)=nkpt
    d(dim2)=nproc/nkpt
  endif    
else if (task.eq.800.or.task.eq.801.or.task.eq.802.or.task.eq.810.or.task.eq.809.or.task.eq.822) then
  i2=nvq
  if (i2.eq.0) i2=nkptnr
  nd=3
  allocate(d(nd))
  d=1
  if (nproc.le.nkptnr) then
    d(dim1)=nproc
  else  
    d(dim1)=nkptnr
    i1=nproc/nkptnr
    if (i1.le.i2) then
      d(dim2)=i1
    else
      d(dim2)=i2
      d(dim3)=nproc/(nkptnr*i2)
    endif
  endif
else
  nd=1
  allocate(d(nd))
  d=nproc
endif  
! overwrite default grid layout
if (lmpigrid) then
  d(1:nd)=mpigrid(1:nd) 
endif
call mpi_grid_initialize(d)
deallocate(d)
return
end