program pp_u
use mod_hdf5
implicit none

character*100 fname
integer nepts,size,ntr_uscrn,it,n,nwann,iw,iwloc,nwloc,i,j
integer n1,n2
integer, allocatable :: vtl_uscrn(:,:)
real(8), allocatable :: w(:)
complex(8), allocatable :: uscrn(:,:,:)
complex(8), allocatable :: ubare(:,:)

complex(8), allocatable :: uscrn_stat(:,:,:)
complex(8), allocatable :: zf(:)
real(8) t1,uavg
character*8 c8
character*3 c3
logical exist
real(8), parameter :: ha2ev = 27.21138386d0
integer vtl(3)
integer nwann_stat
integer, allocatable :: iwann_stat(:)

call hdf5_initialize

open(150,file="pp_u.in",form="FORMATTED",status="OLD")
read(150,*)vtl
read(150,*)nwann_stat
allocate(iwann_stat(nwann_stat))
read(150,*)iwann_stat
close(150)


!nwann=17

fname="uscrn0000.hdf"
call hdf5_read(fname,"/parameters","nepts",nepts)
call hdf5_read(fname,"/parameters","nwann",nwann)
call hdf5_read(fname,"/parameters","size",size)
call hdf5_read(fname,"/parameters","ntr_uscrn",ntr_uscrn)
allocate(vtl_uscrn(3,ntr_uscrn))
call hdf5_read(fname,"/parameters","vtl_uscrn",vtl_uscrn(1,1),(/3,ntr_uscrn/))  

allocate(w(nepts))
allocate(uscrn(nwann,nwann,nepts))
allocate(ubare(nwann,nwann))
do it=1,ntr_uscrn
  if (vtl_uscrn(1,it).eq.vtl(1).and.vtl_uscrn(2,it).eq.vtl(2).and.&
      vtl_uscrn(3,it).eq.vtl(3)) then
    write(c3,'(I3.3)')it
  endif
enddo

do n=0,size-1
  write(fname,'("uscrn",I4.4,".hdf")')n
  inquire(file=trim(fname),exist=exist)
  if (exist) then
    call hdf5_read(fname,"/parameters","nwloc",nwloc)
    do iwloc=1,nwloc
      write(c8,'(I8.8)')iwloc
      call hdf5_read(fname,"/iwloc/"//c8,"iw",iw)
      call hdf5_read(fname,"/iwloc/"//c8,"w",w(iw))
      call hdf5_read(fname,"/iwloc/"//c8//"/"//c3,"uscrn",uscrn(1,1,iw),(/nwann,nwann/))
    enddo
    if (n.eq.0) then
      call hdf5_read(fname,"/iwloc/00000001/"//c3,"ubare",ubare(1,1),(/nwann,nwann/))
    endif
  endif
enddo

allocate(uscrn_stat(nwann_stat,nwann_stat,nepts))
do iw=1,nepts
  do n1=1,nwann_stat
    do n2=1,nwann_stat
      uscrn_stat(n1,n2,iw)=uscrn(iwann_stat(n1),iwann_stat(n2),iw)
    enddo
  enddo
enddo

allocate(zf(2))    
open(150,file="cRPA.dat",status="REPLACE",form="FORMATTED")
write(150,'("# Wannier functions : ",100I4)')iwann_stat
write(150,'("#")')
write(150,'("# Screened U(w=0) matrix")')
write(150,'("#  real part")')
do i=1,nwann_stat
  write(150,'("# ",100F12.6)')(dreal(uscrn_stat(i,j,1)),j=1,nwann_stat)
enddo
write(150,'("#  imag part")')
do i=1,nwann_stat
  write(150,'("# ",100F12.6)')(dimag(uscrn_stat(i,j,1)),j=1,nwann_stat)
enddo
write(150,'("#")')
write(150,'("# Bare U matrix")')
write(150,'("#  real part")')
do i=1,nwann_stat
  write(150,'("# ",100F12.6)')(dreal(ubare(iwann_stat(i),iwann_stat(j))),j=1,nwann_stat)
enddo
write(150,'("#  imag part")')
do i=1,nwann_stat
  write(150,'("# ",100F12.6)')(dimag(ubare(iwann_stat(i),iwann_stat(j))),j=1,nwann_stat)
enddo
write(150,'("#")')
write(150,'("# columns : ")')
write(150,'("#   1 : enery ")')
write(150,'("#   2 : Re(U_avg_diagonal) ")')
write(150,'("#   3 : Im(U_avg_diagonal) ")')
write(150,'("#   4 : Re(U_avg_total) ")')
write(150,'("#   5 : Im(U_avg_total) ")')
write(150,'("#   6 : Re(U_{11}) ")')
write(150,'("#   7 : Im(U_{11}) ")')
write(150,'("#   8 : Re(U_{22}) ")')
write(150,'("#   9 : Im(U_{22}) ")')
write(150,'("#   10 : ... ")')


write(150,'("#")')
do iw=1,nepts
  zf=dcmplx(0.d0,0.d0)
  do n1=1,nwann_stat
    zf(1)=zf(1)+uscrn_stat(n1,n1,iw)/nwann_stat
  enddo
  do n1=1,nwann_stat
    do n2=1,nwann_stat
      zf(2)=zf(2)+uscrn_stat(n1,n2,iw)/nwann_stat/nwann_stat
    enddo
  enddo
  write(150,'(15G18.10)')w(iw)*ha2ev,dreal(zf(1)),dimag(zf(1)),&
    dreal(zf(2)),dimag(zf(2)),(dreal(uscrn_stat(n1,n1,iw)),&
    dimag(uscrn_stat(n1,n1,iw)),n1=1,nwann_stat)
enddo
close(150)

open(150,file="CRPA_U.OUT",status="REPLACE",form="FORMATTED")
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
  uavg=uavg+dreal(uscrn(i,i,1))/nwann
enddo
write(150,'("Average diagonal screened U : ",F12.6)')uavg
uavg=0.d0
do i=1,nwann
  do j=1,nwann
    uavg=uavg+dreal(uscrn(i,j,1))/nwann/nwann
  enddo
enddo
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