subroutine write_u
use modmain
implicit none
integer i,j
real(8) uavg

if (write_chi0_file) then
  call mpi_grid_reduce(uscrnwan(1,1,1),nwann*nwann*nepts,dims=(/dim_k,dim_q/),&
    side=.true.)
  call mpi_grid_reduce(ubarewan(1,1),nwann*nwann,dims=(/dim_k,dim_q/),&
    side=.true.)
else
  call mpi_grid_reduce(uscrnwan(1,1,1),nwann*nwann*nepts,dims=(/dim_q/),&
    side=.true.)
  call mpi_grid_reduce(ubarewan(1,1),nwann*nwann,dims=(/dim_q/),side=.true.)
endif

if (mpi_grid_root()) then
  uscrnwan=ha2ev*uscrnwan/omega/nvq0
  ubarewan=ha2ev*ubarewan/omega/nvq0
  open(150,file='crpa.dat',status='replace',form='formatted')
  do i=1,nepts
    write(150,'(3G18.10)')dreal(lr_w(i))*ha2ev, &
      sum(dreal(uscrnwan(:,:,i)))/nwann/nwann, &
      sum(dimag(uscrnwan(:,:,i)))/nwann/nwann
  enddo
  close(150)
  open(150,file='CRPA_U.OUT',status='replace',form='formatted')
  write(150,'("Screened U matrix")')
  write(150,'("real part")')
  do i=1,nwann
    write(150,'(100F12.6)')(dreal(uscrnwan(i,j,1)),j=1,nwann)
  enddo
  write(150,'("imag part")')
  do i=1,nwann
    write(150,'(100F12.6)')(dimag(uscrnwan(i,j,1)),j=1,nwann)
  enddo
  uavg=0.d0
  do i=1,nwann
    uavg=uavg+dreal(uscrnwan(i,i,1))
  enddo
  uavg=uavg/nwann
  write(150,'("Average diagonal screened U : ",F12.6)')uavg
  uavg=0.d0
  do i=1,nwann
    do j=1,nwann
      uavg=uavg+dreal(uscrnwan(i,j,1))
    enddo
  enddo
  uavg=uavg/nwann/nwann
  write(150,'("Average total screened U : ",F12.6)')uavg 
  write(150,*)
  write(150,'("Bare U matrix")')
  write(150,'("real part")')
  do i=1,nwann
    write(150,'(100F12.6)')(dreal(ubarewan(i,j)),j=1,nwann)
  enddo
  write(150,'("imag part")')
  do i=1,nwann
    write(150,'(100F12.6)')(dimag(ubarewan(i,j)),j=1,nwann)
  enddo
  uavg=0.d0
  do i=1,nwann
    uavg=uavg+dreal(ubarewan(i,i))
  enddo
  uavg=uavg/nwann
  write(150,'("Average diagonal bare U : ",F12.6)')uavg
  uavg=0.d0
  do i=1,nwann
    do j=1,nwann
      uavg=uavg+dreal(ubarewan(i,j))
    enddo
  enddo
  uavg=uavg/nwann/nwann
  write(150,'("Average total bare U : ",F12.6)')uavg
  close(150)  
endif
  



return
end