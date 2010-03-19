subroutine write_u
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
call mpi_grid_reduce(ubarewan(1,1),nwann*nwann,dims=(/dim_q/),side=.true.)
ubarewan=ha2ev*ubarewan/omega/nvq0

fuscrn="uscrn.hdf5"
if (mpi_grid_root()) then
  call hdf5_create_file(trim(fuscrn))
  !call hdf5_create_group(trim(fuscrn),"/","parameters")
  !call mdf5_write(fuscrn,"/parameters","ntr",ntr_uscrn)
  !call mdf5_write(fuscrn,"/parameters","vtl",vtl_uscrn,(/3,ntr_uscrn/))
  !call mdf5_write(fuscrn,"/parameters","ivtit",ivtit_uscrn,(/3,ntr_uscrn/)) 
  call hdf5_create_group(trim(fuscrn),"/","iw")
  do i=1,nepts
    write(c8,'(I8.8)')i
    call hdf5_create_group(trim(fuscrn),"/iw",c8)
    do it=1,ntr_uscrn
      write(c3,'(I3.3)')it
      call hdf5_create_group(trim(fuscrn),"/iw"//"/"//c8,c3)
    enddo
  enddo
endif
if (mpi_grid_side(dims=(/dim_k/))) then
  do i=0,mpi_grid_size(dim_k)-1
    if (mpi_grid_x(dim_k).eq.i) then
      do iwloc=1,nwloc
        iw=mpi_grid_map(nepts,dim_k,loc=iwloc)
        write(c8,'(I8.8)')iw
        do it=1,ntr_uscrn
          write(c3,'(I3.3)')it
          call hdf5_write(fuscrn,"/iw"//"/"//c8//"/"//c3,"uscrn",uscrnwan(1,1,it,iwloc),&
            (/nwann,nwann/))
        enddo
      enddo
    endif
    call mpi_grid_barrier(dims=(/dim_k/))
  enddo
endif

        

!if (mpi_grid_root()) then
!  uscrnwan=ha2ev*uscrnwan/omega/nvq0
!  ubarewan=ha2ev*ubarewan/omega/nvq0
!!  open(150,file='crpa.dat',status='replace',form='formatted')
!!  do i=1,nepts
!!    write(150,'(3G18.10)')dreal(lr_w(i))*ha2ev, &
!!      sum(dreal(uscrnwan(:,:,i)))/nwann/nwann, &
!!      sum(dimag(uscrnwan(:,:,i)))/nwann/nwann
!!  enddo
!!  close(150)
!  open(150,file='CRPA_U.OUT',status='replace',form='formatted')
!  write(150,'("Screened U matrix")')
!  write(150,'("real part")')
!  do i=1,nwann
!    write(150,'(100F12.6)')(dreal(uscrnwan(i,j,ivtit_uscrn(0,0,0),1)),j=1,nwann)
!  enddo
!  write(150,'("imag part")')
!  do i=1,nwann
!    write(150,'(100F12.6)')(dimag(uscrnwan(i,j,ivtit_uscrn(0,0,0),1)),j=1,nwann)
!  enddo
!  uavg=0.d0
!  do i=1,nwann
!    uavg=uavg+dreal(uscrnwan(i,i,ivtit_uscrn(0,0,0),1))
!  enddo
!  uavg=uavg/nwann
!  write(150,'("Average diagonal screened U : ",F12.6)')uavg
!  uavg=0.d0
!  do i=1,nwann
!    do j=1,nwann
!      uavg=uavg+dreal(uscrnwan(i,j,ivtit_uscrn(0,0,0),1))
!    enddo
!  enddo
!  uavg=uavg/nwann/nwann
!  write(150,'("Average total screened U : ",F12.6)')uavg 
!  write(150,*)
!  write(150,'("Bare U matrix")')
!  write(150,'("real part")')
!  do i=1,nwann
!    write(150,'(100F12.6)')(dreal(ubarewan(i,j)),j=1,nwann)
!  enddo
!  write(150,'("imag part")')
!  do i=1,nwann
!    write(150,'(100F12.6)')(dimag(ubarewan(i,j)),j=1,nwann)
!  enddo
!  uavg=0.d0
!  do i=1,nwann
!    uavg=uavg+dreal(ubarewan(i,i))
!  enddo
!  uavg=uavg/nwann
!  write(150,'("Average diagonal bare U : ",F12.6)')uavg
!  uavg=0.d0
!  do i=1,nwann
!    do j=1,nwann
!      uavg=uavg+dreal(ubarewan(i,j))
!    enddo
!  enddo
!  uavg=uavg/nwann/nwann
!  write(150,'("Average total bare U : ",F12.6)')uavg
!  close(150)  
!endif
  



return
end