subroutine write_ubare
use modmain
implicit none
character*100 fubare
character*8 c8
character*3 c3
integer it

write(*,*)'ubare=',ubarewan(1,1,14)

call mpi_grid_reduce(ubarewan(1,1,1),nwann*nwann*ntr_uscrn,dims=(/dim_b,dim_q/),&
  side=.true.)
ubarewan=ha2ev*ubarewan !/omega/nkptnr

if (mpi_grid_side(dims=(/dim_k/)).and.mpi_grid_x(dim_k).eq.0) then
  fubare="ubare.hdf"
  call hdf5_create_file(trim(fubare))
  call hdf5_create_group(trim(fubare),"/","parameters")
  call hdf5_create_group(trim(fubare),"/","it")
  call hdf5_write(fubare,"/parameters","nwann",nwann)
  call hdf5_write(fubare,"/parameters","ntr_uscrn",ntr_uscrn)
  call hdf5_write(fubare,"/parameters","vtl_uscrn",vtl_uscrn(1,1),(/3,ntr_uscrn/))  
  if (mpi_grid_x(dim_k).eq.0) then
    do it=1,ntr_uscrn
      write(c3,'(I3.3)')it
      call hdf5_create_group(trim(fubare),"/it",c3)
      call hdf5_write(fubare,"/it/"//c3,"ubare",ubarewan(1,1,it),(/nwann,nwann/))
    enddo  
  endif
endif

return
end