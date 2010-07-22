program pp_u
use mod_hdf5
implicit none

character*100 fname
integer nw,size,ntr_uscrn,it,n,nwann,iw,iwloc,nwloc,i,j
integer n1,n2
real(8), allocatable :: w(:)
complex(8), allocatable :: uscrn(:)

complex(8), allocatable :: uscrn_(:,:,:)
complex(8) zf(2)
real(8) t1,uavg
character*8 c8
logical exist
real(8), parameter :: ha2ev = 27.21138386d0
integer vtl_(3)
integer nwann_
integer, allocatable :: iwann_(:)
integer nmegqwan
integer, allocatable :: imegqwan(:,:)
integer, allocatable :: iwm(:,:)
call hdf5_initialize

open(150,file="pp_u.in",form="FORMATTED",status="OLD")
read(150,*)vtl_
read(150,*)nwann_
allocate(iwann_(nwann_))
read(150,*)iwann_
close(150)

allocate(iwm(nwann_,nwann_))



fname="uscrn0000.hdf"
inquire(file=trim(fname),exist=exist)
if (exist) then
  call hdf5_read(fname,"/parameters","nw",nw)
  call hdf5_read(fname,"/parameters","nwann",nwann)
  call hdf5_read(fname,"/parameters","size",size)
  call hdf5_read(fname,"/parameters","nmegqwan",nmegqwan)
  allocate(imegqwan(5,nmegqwan))
  call hdf5_read(fname,"/parameters","imegqwan",imegqwan(1,1),(/5,nmegqwan/))  
  allocate(w(nw))
  allocate(uscrn(nmegqwan))
  allocate(uscrn_(nwann_,nwann_,nw))
  uscrn_=dcmplx(0.d0,0.d0)
! create i -> {nn'T} mapping
  do i=1,nmegqwan
    do n1=1,nwann_
      do n2=1,nwann_
        if (iwann_(n1).eq.imegqwan(1,i).and.&
            iwann_(n2).eq.imegqwan(2,i).and.&
            all(vtl_(:).eq.imegqwan(3:5,i))) then
            iwm(n1,n2)=i
        endif
      enddo
    enddo
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
        call hdf5_read(fname,"/iwloc/"//c8,"uscrn",uscrn(1),(/nmegqwan/))
        do n1=1,nwann_
          do n2=1,nwann_
            uscrn_(n1,n2,iw)=uscrn(iwm(n1,n2))*ha2ev
          enddo
        enddo
      enddo
    endif
  enddo
  
  open(150,file="cRPA.dat",status="REPLACE",form="FORMATTED")
  write(150,'("# Wannier functions : ",100I4)')iwann_
  write(150,'("#")')
  write(150,'("# Screened U(w=0) matrix")')
  write(150,'("#  real part")')
  do i=1,nwann_
    write(150,'("# ",100F12.6)')(dreal(uscrn_(i,j,1)),j=1,nwann_)
  enddo
  write(150,'("#  imag part")')
  do i=1,nwann_
    write(150,'("# ",100F12.6)')(dimag(uscrn_(i,j,1)),j=1,nwann_)
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
  do iw=1,nw
    zf=dcmplx(0.d0,0.d0)
    do n1=1,nwann_
      zf(1)=zf(1)+uscrn_(n1,n1,iw)/nwann_
    enddo
    do n1=1,nwann_
      do n2=1,nwann_
        zf(2)=zf(2)+uscrn_(n1,n2,iw)/nwann_/nwann_
      enddo
    enddo
    write(150,'(100G18.10)')w(iw)*ha2ev,dreal(zf(1)),dimag(zf(1)),&
      dreal(zf(2)),dimag(zf(2)),(dreal(uscrn_(n1,n1,iw)),&
      dimag(uscrn_(n1,n1,iw)),n1=1,nwann_)
  enddo
  close(150)
  
!  open(150,file="CRPA_U.OUT",status="REPLACE",form="FORMATTED")
!  write(150,'("Screened U matrix")')
!  write(150,'("real part")')
!  do i=1,nwann
!    write(150,'(100F12.6)')(dreal(uscrn(i,j,1)),j=1,nwann)
!  enddo
!  write(150,'("imag part")')
!  do i=1,nwann
!    write(150,'(100F12.6)')(dimag(uscrn(i,j,1)),j=1,nwann)
!  enddo
!  uavg=0.d0
!  do i=1,nwann
!    uavg=uavg+dreal(uscrn(i,i,1))/nwann
!  enddo
!  write(150,'("Average diagonal screened U : ",F12.6)')uavg
!  uavg=0.d0
!  do i=1,nwann
!    do j=1,nwann
!      uavg=uavg+dreal(uscrn(i,j,1))/nwann/nwann
!    enddo
!  enddo
!  write(150,'("Average total screened U : ",F12.6)')uavg 
!  close(150)
endif
call hdf5_finalize
return
end

subroutine pstop
stop
end