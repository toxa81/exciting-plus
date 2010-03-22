program pp_u
use mod_hdf5
implicit none

character*100 fname
integer nepts,size,ntr_uscrn,it,n,nwann,iw,iwloc,nwloc,i,j
integer, allocatable :: vtl_uscrn(:,:)
real(8), allocatable :: w(:)
complex(8), allocatable :: uscrn(:,:,:)
real(8) t1,uavg
character*8 c8
character*3 c3
real(8), parameter :: ha2ev = 27.21138386d0

call hdf5_initialize

nwann=8

fname="uscrn0000.hdf"
call hdf5_read(fname,"/parameters","nepts",nepts)
!call hdf5_read(fname,"/parameters","nwloc",nwloc)
!call hdf5_read(fname,"/parameters","x",mpi_grid_x(dim_k))
call hdf5_read(fname,"/parameters","size",size)
call hdf5_read(fname,"/parameters","ntr_uscrn",ntr_uscrn)
allocate(vtl_uscrn(3,ntr_uscrn))
call hdf5_read(fname,"/parameters","vtl_uscrn",vtl_uscrn(1,1),(/3,ntr_uscrn/))  

allocate(w(nepts))
allocate(uscrn(nwann,nwann,nepts))
do it=1,ntr_uscrn
  if (vtl_uscrn(1,it).eq.0.and.vtl_uscrn(2,it).eq.0.and.vtl_uscrn(3,it).eq.0) then
    write(c3,'(I3.3)')it
  endif
enddo
do n=0,size-1
  write(fname,'("uscrn",I4.4,".hdf")')n
  call hdf5_read(fname,"/parameters","nwloc",nwloc)
  do iwloc=1,nwloc
    write(c8,'(I8.8)')iwloc
    call hdf5_read(fname,"/iwloc/"//c8,"iw",iw)
    call hdf5_read(fname,"/iwloc/"//c8,"w",w(iw))
    call hdf5_read(fname,"/iwloc/"//c8//"/"//c3,"uscrn",uscrn(1,1,iw))
  enddo
enddo
    
    
open(150,file='crpa.dat',status='replace',form='formatted')
do i=1,nepts
  write(150,'(3G18.10)')w(i)*ha2ev, &
    sum(dreal(uscrn(:,:,i)))/nwann/nwann, &
    sum(dimag(uscrn(:,:,i)))/nwann/nwann
enddo
close(150)

open(150,file='CRPA_U.OUT',status='replace',form='formatted')
write(150,'("Screened U matrix")')
write(150,'("real part")')
do i=1,nwann
  write(150,'(100F12.6)')(dreal(uscrn(i,j,1)),j=1,nwann)
enddo
write(150,'("imag part")')
do i=1,nwann
  write(150,'(100F12.6)')(dimag(uscrn(i,j,1)),j=1,nwann)
enddo
uavg=0.d0
do i=1,nwann
  uavg=uavg+dreal(uscrn(i,i,1))
enddo
uavg=uavg/nwann
write(150,'("Average diagonal screened U : ",F12.6)')uavg
uavg=0.d0
do i=1,nwann
  do j=1,nwann
    uavg=uavg+dreal(uscrn(i,j,1))
  enddo
enddo
uavg=uavg/nwann/nwann
write(150,'("Average total screened U : ",F12.6)')uavg 
!write(150,*)
!write(150,'("Bare U matrix")')
!write(150,'("real part")')
!do i=1,nwann
!  write(150,'(100F12.6)')(dreal(ubarewan(i,j)),j=1,nwann)
!enddo
!write(150,'("imag part")')
!do i=1,nwann
!  write(150,'(100F12.6)')(dimag(ubarewan(i,j)),j=1,nwann)
!enddo
!uavg=0.d0
!do i=1,nwann
!  uavg=uavg+dreal(ubarewan(i,i))
!enddo
!uavg=uavg/nwann
!write(150,'("Average diagonal bare U : ",F12.6)')uavg
!uavg=0.d0
!do i=1,nwann
!  do j=1,nwann
!    uavg=uavg+dreal(ubarewan(i,j))
!  enddo
!enddo
!uavg=uavg/nwann/nwann
!write(150,'("Average total bare U : ",F12.6)')uavg
!close(150)  
!
!



call hdf5_finalize
return
end