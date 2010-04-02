subroutine write_uscrn
use modmain
implicit none
integer i,j,nwloc,it,iwloc,iw
real(8) uavg
character*100 fuscrn
character*8 c8
character*3 c3

nwloc=mpi_grid_map(nepts,dim_k)

call mpi_grid_reduce(uscrnwan(1,1,1,1),nwann*nwann*ntr_uscrn*nwloc,&
  dims=(/dim_b,dim_q/))
uscrnwan=ha2ev*uscrnwan/omega/nvq0

if (mpi_grid_side(dims=(/dim_k/)).and.nwloc.gt.0) then
  write(fuscrn,'("uscrn",I4.4,".hdf")')mpi_grid_x(dim_k)
  call hdf5_create_file(trim(fuscrn))
  call hdf5_create_group(trim(fuscrn),"/","iwloc")
  call hdf5_create_group(trim(fuscrn),"/","parameters")
  call hdf5_write(fuscrn,"/parameters","nwann",nwann)
  call hdf5_write(fuscrn,"/parameters","nepts",nepts)
  call hdf5_write(fuscrn,"/parameters","nwloc",nwloc)
  call hdf5_write(fuscrn,"/parameters","x",mpi_grid_x(dim_k))
  call hdf5_write(fuscrn,"/parameters","size",mpi_grid_size(dim_k))
  call hdf5_write(fuscrn,"/parameters","ntr_uscrn",ntr_uscrn)
  call hdf5_write(fuscrn,"/parameters","vtl_uscrn",vtl_uscrn(1,1),(/3,ntr_uscrn/))  
  do iwloc=1,nwloc
    iw=mpi_grid_map(nepts,dim_k,loc=iwloc)
    write(c8,'(I8.8)')iwloc
    call hdf5_create_group(trim(fuscrn),"/iwloc",c8)
    call hdf5_write(fuscrn,"/iwloc/"//c8,"iw",iw)
    call hdf5_write(fuscrn,"/iwloc/"//c8,"w",dreal(lr_w(iw)))
    do it=1,ntr_uscrn
      write(c3,'(I3.3)')it
      call hdf5_create_group(trim(fuscrn),"/iwloc/"//c8,c3)
      call hdf5_write(fuscrn,"/iwloc/"//c8//"/"//c3,"uscrn",&
        uscrnwan(1,1,it,iwloc),(/nwann,nwann/))
    enddo
  enddo
endif
return
end