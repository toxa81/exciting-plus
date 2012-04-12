subroutine initmpigrid
use modmain
use mod_addons_q
implicit none
integer, allocatable :: d(:)
integer i1,i2,nd,i
select case(task)
  case(0,1,20,21,805,822,700,701,702)
    nd=2
    allocate(d(nd)); d(:)=1
    if (nproc.le.nkpt) then
      d(dim1)=nproc
    else
      d(dim1)=nkpt
      d(dim2)=nproc/nkpt
    endif
! linear response runs on 3D grid
  case(800,801)
    nd=3
    allocate(d(nd)); d(:)=1
    i2=nvq
    if (i2.eq.0) i2=nkptnr
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
  case default
    nd=1
    allocate(d(nd))
    d=nproc
end select  
! overwrite default grid layout
if (lmpigrid) then
  deallocate(d)
  allocate(d(mpigrid_ndim))
  d(1:mpigrid_ndim)=mpigrid(1:mpigrid_ndim)
endif
call mpi_grid_initialize(d)
if (mpi_grid_root()) then
  write(*,*)
  write(*,'("[initmpigrid] mpi grid size : ",10I8)')(mpi_grid_dim_size(i),i=1,mpi_grid_nd)
endif
deallocate(d)
return
end
