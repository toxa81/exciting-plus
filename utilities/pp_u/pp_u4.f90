program pp_u
use mod_hdf5
implicit none
integer i,j,is1,is2,n1,n2,n,it,iw,nwloc,iwloc
character*100 str,fname,fout
integer mode,vtl_(3),nw,nwann,size,ngq,ngridk(3)
integer nwann_
integer, allocatable :: iwann_(:)
integer nmegqwan
integer, allocatable :: imegqwan(:,:)
integer ntmegqwan
integer, allocatable :: itmegqwan(:,:)
integer tlim1(2,3)
integer, allocatable :: idxt1(:,:,:,:,:) 
integer tlim2(2,3)
integer, allocatable :: idxt2(:,:,:) 
real(8), allocatable :: w(:)
complex(8), allocatable :: u4(:,:,:)
complex(8), allocatable :: uscrn(:,:,:)
complex(8), allocatable :: jscrn(:,:,:)
complex(8), allocatable :: ujmtrx(:,:)
character*8 c1,c2,c3,c4
complex(8) uav,jav,zeps
logical exist
real(8), parameter :: ha2ev = 27.21138386d0

call hdf5_initialize
open(150,file="pp_u.in",form="FORMATTED",status="OLD")
read(150,'(A)')str
if (trim(adjustl(str)).eq."u2") then
  mode=0
  read(150,*)nwann_
  allocate(iwann_(nwann_))
  read(150,*)iwann_
  read(150,*)vtl_
  close(150)
!else if (trim(adjustl(str)).eq."offsite") then
!  mode=1
!else if (trim(adjustl(str)).eq."list") then
!  mode=2
!  read(150,*)nlist_
!  allocate(ilist_(nlist_))
!  read(150,*)ilist_
endif

fname="uscrn0000.hdf5"
inquire(file=trim(fname),exist=exist)
if (.not.exist) then
  write(*,'("File uscrn0000.hdf5 not found.")')
  return
endif
call hdf5_read(fname,"/parameters","nw",nw)
call hdf5_read(fname,"/parameters","nwann",nwann)
call hdf5_read(fname,"/parameters","size",size)
call hdf5_read(fname,"/parameters","nmegqwan",nmegqwan)
allocate(imegqwan(5,nmegqwan))
call hdf5_read(fname,"/parameters","imegqwan",imegqwan(1,1),(/5,nmegqwan/))  
call hdf5_read(fname,"/parameters","ngq",ngq)
call hdf5_read(fname,"/parameters","ngridk",ngridk(1),(/3/))  
call hdf5_read(fname,"/parameters","ntmegqwan",ntmegqwan)
allocate(itmegqwan(3,ntmegqwan))
call hdf5_read(fname,"/parameters","itmegqwan",itmegqwan(1,1),(/3,ntmegqwan/))
! find translation limits and build reverce ((tx,ty,tz) -> it) index
tlim1(1,1)=minval(imegqwan(3,:))
tlim1(2,1)=maxval(imegqwan(3,:))
tlim1(1,2)=minval(imegqwan(4,:))
tlim1(2,2)=maxval(imegqwan(4,:))
tlim1(1,3)=minval(imegqwan(5,:))
tlim1(2,3)=maxval(imegqwan(5,:))
allocate(idxt1(nwann,nwann,tlim1(1,1):tlim1(2,1),tlim1(1,2):tlim1(2,2),&
  tlim1(1,3):tlim1(2,3)))
idxt1=-100
do i=1,nmegqwan
  idxt1(imegqwan(1,i),imegqwan(2,i),imegqwan(3,i),imegqwan(4,i),&
    imegqwan(5,i))=i
enddo
tlim2(1,1)=minval(itmegqwan(1,:))
tlim2(2,1)=maxval(itmegqwan(1,:))
tlim2(1,2)=minval(itmegqwan(2,:))
tlim2(2,2)=maxval(itmegqwan(2,:))
tlim2(1,3)=minval(itmegqwan(3,:))
tlim2(2,3)=maxval(itmegqwan(3,:))
allocate(idxt2(tlim2(1,1):tlim2(2,1),tlim2(1,2):tlim2(2,2),tlim2(1,3):tlim2(2,3)))
idxt2=-100
do i=1,ntmegqwan
  idxt2(itmegqwan(1,i),itmegqwan(2,i),itmegqwan(3,i))=i
enddo
allocate(w(nw))
allocate(u4(nmegqwan,nmegqwan,ntmegqwan))
write(c1,'(I6)')ngridk(1)
write(c2,'(I6)')ngridk(2)
write(c3,'(I6)')ngridk(3)
write(c4,'(I8)')ngq
fout="cRPA__"//trim(adjustl(c1))//"x"//trim(adjustl(c2))//"x"//&
  trim(adjustl(c3))//"_k__"//trim(adjustl(c4))//"_G__.dat"  
zeps=dcmplx(0.d0,1d-12)
if (mode.eq.0) then
! check that we have the required translation
  if (idxt2(vtl_(1),vtl_(2),vtl_(3)).lt.0) then
    write(*,'("Error: translation ",3I4," is not in list")')vtl_
  endif
  allocate(uscrn(nwann_,nwann_,nw))
  allocate(jscrn(nwann_,nwann_,nw))
  allocate(ujmtrx(2*nwann_,2*nwann_))
  uscrn=dcmplx(0.d0,0.d0)
  jscrn=dcmplx(0.d0,0.d0)
  do n=0,size-1
    write(fname,'("uscrn",I4.4,".hdf5")')n
    inquire(file=trim(fname),exist=exist)
    if (exist) then
      call hdf5_read(fname,"/parameters","nwloc",nwloc)
      do iwloc=1,nwloc
        write(c1,'(I8.8)')iwloc
        call hdf5_read(fname,"/iwloc/"//c1,"iw",iw)
        call hdf5_read(fname,"/iwloc/"//c1,"w",w(iw))
        do it=1,ntmegqwan
          write(c2,'("t",I7.7)')it
          call hdf5_read(fname,"/iwloc/"//c1//"/"//c2,"u4",u4(1,1,it),&
            (/nmegqwan,nmegqwan/))
        enddo
! take 2-index U
        do n1=1,nwann_
          do n2=1,nwann_
            uscrn(n1,n2,iw)=&
              u4(idxt1(iwann_(n1),iwann_(n1),0,0,0),&
                 idxt1(iwann_(n2),iwann_(n2),0,0,0),&
                 idxt2(vtl_(1),vtl_(2),vtl_(3)))*ha2ev
            jscrn(n1,n2,iw)=&
              u4(idxt1(iwann_(n1),iwann_(n2),0,0,0),&
                 idxt1(iwann_(n2),iwann_(n1),0,0,0),&
                 idxt2(vtl_(1),vtl_(2),vtl_(3)))*ha2ev
          enddo !n2
        enddo !n1
      enddo !iwloc
    endif !exist
  enddo !n 
! make DMFT "density-density" input matrix    
  do is1=1,2
    do is2=1,2
      do i=1,nwann_
        do j=1,nwann_
          if (is1.eq.is2) then
            ujmtrx((is1-1)*nwann_+i,(is2-1)*nwann_+j)=uscrn(i,j,1)-jscrn(i,j,1)
          else
            ujmtrx((is1-1)*nwann_+i,(is2-1)*nwann_+j)=uscrn(i,j,1)
          endif
        enddo
      enddo
    enddo
  enddo
! write to file   
  open(150,file=trim(adjustl(fout)),status="REPLACE",form="FORMATTED")
  write(150,'("# Wannier functions : ",100I4)')iwann_
  write(150,'("#")')
  write(150,'("# Screened U(w=0) matrix")')
  write(150,'("# real part ")')  
  do i=1,nwann_
    write(150,'("# ",100F16.8)')(dreal(uscrn(i,j,1)),j=1,nwann_)
  enddo
  write(150,'("# image part ")')  
  do i=1,nwann_
    write(150,'("# ",100F16.8)')(dimag(uscrn(i,j,1)-zeps),j=1,nwann_)
  enddo
  write(150,'("#")')    
  write(150,'("# Screened J(w=0) matrix")')
  write(150,'("# real part ")')  
  do i=1,nwann_
    write(150,'("# ",100F16.8)')(dreal(jscrn(i,j,1)),j=1,nwann_)
  enddo
  write(150,'("# image part ")')  
  do i=1,nwann_
    write(150,'("# ",100F16.8)')(dimag(jscrn(i,j,1)-zeps),j=1,nwann_)
  enddo
  write(150,'("#")')    
  write(150,'("# UJ matrix")')
  write(150,'("# real part ")')    
  do i=1,2*nwann_
    write(150,'("# ",100F16.8)')(dreal(ujmtrx(i,j)),j=1,2*nwann_)
  enddo
  write(150,'("# image part ")')    
  do i=1,2*nwann_
    write(150,'("# ",100F16.8)')(dimag(ujmtrx(i,j)-zeps),j=1,2*nwann_)
  enddo
  write(150,'("#")')
  write(150,'("# columns (units are eV) : ")')
  write(150,'("#   1  : energy ")')
  write(150,'("#   2  : Re (average U) ")')
  write(150,'("#   3  : Im (average U) ")')
  write(150,'("#   4  : Re (average J) ")')
  write(150,'("#   5  : Im (average J) ")')
  write(150,'("#   6  : Re (U_{11}) ")')
  write(150,'("#   7  : Im (U_{11}) ")')
  write(150,'("#   8  : Re (J_{11}) ")')
  write(150,'("#   9  : Im (J_{11}) ")')
  write(150,'("#   10 : Re (U_{22}) ")')
  write(150,'("#   11 : ... ")')
  
  write(150,'("#")')
  do iw=1,nw
    uav=dcmplx(0.d0,0.d0)
    do n1=1,nwann_
      do n2=1,nwann_
        uav=uav+uscrn(n1,n2,iw)
      enddo
    enddo
    uav=uav/nwann_/nwann_
    jav=dcmplx(0.d0,0.d0)
    do n1=1,nwann_
      do n2=1,nwann_
        if (n1.ne.n2) jav=jav+jscrn(n1,n2,iw)
      enddo
    enddo
    if (nwann_.gt.1) jav=jav/nwann_/(nwann_-1)
    write(150,'(100F16.8)')w(iw)*ha2ev,dreal(uav),dimag(uav-zeps),&
      dreal(jav),dimag(jav-zeps),(dreal(uscrn(n,n,iw)),dimag(uscrn(n,n,iw)-zeps),&
      dreal(jscrn(n,n,iw)),dimag(jscrn(n,n,iw)-zeps),n=1,nwann_)
  enddo
  close(150)
endif !mode.eq.0

!  if (mode.eq.2) then
!    allocate(uscrn2_(nlist_,nw))
!    uscrn2_=dcmplx(0.d0,0.d0)
!    do n=0,size-1
!      write(fname,'("uscrn",I4.4,".hdf")')n
!      inquire(file=trim(fname),exist=exist)
!      if (exist) then
!        call hdf5_read(fname,"/parameters","nwloc",nwloc)
!        do iwloc=1,nwloc
!          write(c8,'(I8.8)')iwloc
!          call hdf5_read(fname,"/iwloc/"//c8,"iw",iw)
!          call hdf5_read(fname,"/iwloc/"//c8,"w",w(iw))
!          call hdf5_read(fname,"/iwloc/"//c8,"uscrn",uscrn(1),(/nmegqwan/))
!          do i=1,nlist_
!            uscrn2_(i,iw)=uscrn(ilist_(i))*ha2ev
!          enddo
!        enddo
!      endif
!    enddo
!    
!    open(150,file=trim(adjustl(fout)),status="REPLACE",form="FORMATTED")
!    write(150,'("#")')
!    write(150,'("# columns : ")')
!    write(150,'("#   1 : energy ")')
!    write(150,'("#   2 : Re(U(first transition in list)) ")')
!    write(150,'("#   3 : Im(U(first transition in list)) ")')
!    write(150,'("#   4 : Re(U(second transition in list)) ")')
!    write(150,'("#   5 : Im(U(second transition in list)) ")')    
!    write(150,'("#   5 : ... ")')
!    write(150,'("#")')
!    do iw=1,nw
!      write(150,'(100G18.10)')w(iw)*ha2ev,(dreal(uscrn2_(i,iw)),&
!        dimag(uscrn2_(i,iw)),i=1,nlist_)
!    enddo
!    close(150)
!  endif
call hdf5_finalize
end

subroutine pstop
stop
end