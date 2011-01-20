program pp_u
use mod_hdf5
implicit none
integer i,j,is1,is2,n1,n2,n,it,iw,nwloc,iwloc
integer i1,i2,j1,j2
character*100 str,fname,fout
integer mode,vtl_(3),nw,nwantot,size,ngq,ngridk(3)
integer nwann_
integer, allocatable :: iwann_(:)
integer nu4
integer, allocatable :: iu4(:,:) 
integer nwt
integer, allocatable :: iwt(:,:)
integer ntr
integer, allocatable :: vtr(:,:)
integer tlim1(2,3)
integer, allocatable :: iwtidx(:,:,:,:,:) 
integer tlim2(2,3)
integer, allocatable :: vtridx(:,:,:) 
real(8), allocatable :: w(:)
complex(8), allocatable :: u4mtrx(:,:,:)
complex(8), allocatable :: uscrn(:,:,:)
complex(8), allocatable :: jscrn(:,:,:)
complex(8), allocatable :: ujmtrx(:,:)
complex(8), allocatable :: u4(:,:)
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
else if (trim(adjustl(str)).eq."u4") then
  mode=1
  read(150,*)nu4
  allocate(iu4(13,nu4))
  do i=1,nu4
    read(150,*)iu4(:,i)
  enddo
!else if (trim(adjustl(str)).eq."list") then
!  mode=2
!  read(150,*)nlist_
!  allocate(ilist_(nlist_))
!  read(150,*)ilist_
endif





fname="u4_0000.hdf5"
inquire(file=trim(fname),exist=exist)
if (.not.exist) then
  write(*,'("File u4_0000.hdf5 not found.")')
  return
endif
call hdf5_read(fname,"/parameters","nw",nw)
call hdf5_read(fname,"/parameters","nwantot",nwantot)
call hdf5_read(fname,"/parameters","size",size)
call hdf5_read(fname,"/parameters","nwt",nwt)
allocate(iwt(5,nwt))
call hdf5_read(fname,"/parameters","iwt",iwt(1,1),(/5,nwt/))  
call hdf5_read(fname,"/parameters","ngq",ngq)
call hdf5_read(fname,"/parameters","ngridk",ngridk(1),(/3/))  
call hdf5_read(fname,"/parameters","ntr",ntr)
allocate(vtr(3,ntr))
call hdf5_read(fname,"/parameters","vtr",vtr(1,1),(/3,ntr/))
! find translation limits and build reverce ((tx,ty,tz) -> it) index
do i=1,3
  tlim1(1,i)=minval(iwt(2+i,:))
  tlim1(2,i)=maxval(iwt(2+i,:))
enddo
allocate(iwtidx(nwantot,nwantot,tlim1(1,1):tlim1(2,1),tlim1(1,2):tlim1(2,2),&
  tlim1(1,3):tlim1(2,3)))
iwtidx=-1
do i=1,nwt
  iwtidx(iwt(1,i),iwt(2,i),iwt(3,i),iwt(4,i),iwt(5,i))=i
enddo
do i=1,3
  tlim2(1,i)=minval(vtr(i,:))
  tlim2(2,i)=maxval(vtr(i,:))
enddo
allocate(vtridx(tlim2(1,1):tlim2(2,1),tlim2(1,2):tlim2(2,2),tlim2(1,3):tlim2(2,3)))
vtridx=-1
do i=1,ntr
  vtridx(vtr(1,i),vtr(2,i),vtr(3,i))=i
enddo
allocate(w(nw))
allocate(u4mtrx(nwt,nwt,ntr))
write(c1,'(I6)')ngridk(1)
write(c2,'(I6)')ngridk(2)
write(c3,'(I6)')ngridk(3)
write(c4,'(I8)')ngq
fout="cRPA__"//trim(adjustl(c1))//"x"//trim(adjustl(c2))//"x"//&
  trim(adjustl(c3))//"_k__"//trim(adjustl(c4))//"_G__.dat"  
zeps=dcmplx(0.d0,1d-12)
if (mode.eq.0) then
! check that we have the required translation
  if (vtridx(vtl_(1),vtl_(2),vtl_(3)).lt.0) then
    write(*,'("Error: translation ",3I4," is not in list")')vtl_
  endif
  allocate(uscrn(nwann_,nwann_,nw))
  allocate(jscrn(nwann_,nwann_,nw))
  allocate(ujmtrx(2*nwann_,2*nwann_))
  uscrn=dcmplx(0.d0,0.d0)
  jscrn=dcmplx(0.d0,0.d0)
  do n=0,size-1
    write(fname,'("u4_",I4.4,".hdf5")')n
    inquire(file=trim(fname),exist=exist)
    if (exist) then
      call hdf5_read(fname,"/parameters","nwloc",nwloc)
      do iwloc=1,nwloc
        write(c1,'(I8.8)')iwloc
        call hdf5_read(fname,"/iwloc/"//c1,"iw",iw)
        call hdf5_read(fname,"/iwloc/"//c1,"w",w(iw))
        do it=1,ntr
          write(c2,'("t",I7.7)')it
          call hdf5_read(fname,"/iwloc/"//c1//"/"//c2,"u4",u4mtrx(1,1,it),&
            (/nwt,nwt/))
        enddo
! take 2-index U
        do n1=1,nwann_
          do n2=1,nwann_
            i1=iwtidx(iwann_(n1),iwann_(n1),0,0,0)
            i2=iwtidx(iwann_(n2),iwann_(n2),0,0,0)
            j1=iwtidx(iwann_(n1),iwann_(n2),0,0,0)
            j2=iwtidx(iwann_(n2),iwann_(n1),0,0,0)
            it=vtridx(vtl_(1),vtl_(2),vtl_(3))
            uscrn(n1,n2,iw)=u4mtrx(i1,i2,it)*ha2ev
            jscrn(n1,n2,iw)=u4mtrx(j1,j2,it)*ha2ev
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
  write(150,'("# U(w=0) matrix")')
  write(150,'("# real part ")')  
  do i=1,nwann_
    write(150,'("# ",100F16.8)')(dreal(uscrn(i,j,1)),j=1,nwann_)
  enddo
  write(150,'("# image part ")')  
  do i=1,nwann_
    write(150,'("# ",100F16.8)')(dimag(uscrn(i,j,1)-zeps),j=1,nwann_)
  enddo
  write(150,'("#")')    
  write(150,'("# J(w=0) matrix")')
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
  write(150,'("#   1 : energy ")')
  write(150,'("#   2 : Re (average U) ")')
  write(150,'("#   3 : Im (average U) ")')
  write(150,'("#   4 : Re (average J) ")')
  write(150,'("#   5 : Im (average J) ")')
  write(150,'("#   6 : Re (U_{11}) ")')
  write(150,'("#   7 : Im (U_{11}) ")')
  write(150,'("#   8 : Re (U_{22}) ")')
  write(150,'("#   9 : Im (U_{22}) ")')
  write(150,'("#  10 : ... ")')
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
      dreal(jav),dimag(jav-zeps),(dreal(uscrn(n,n,iw)),&
        dimag(uscrn(n,n,iw)-zeps),n=1,nwann_)
  enddo
  close(150)
endif !mode.eq.0

if (mode.eq.1) then
! check that we have the required translation
  do i=1,nu4
    if (vtridx(iu4(11,i),iu4(12,i),iu4(13,i)).lt.0) then
      write(*,'("Error: translation ",3I4," is not in list")')iu4(11:13,i)
    endif
  enddo
  allocate(u4(nu4,nw))
  u4=dcmplx(0.d0,0.d0)
  do n=0,size-1
    write(fname,'("uscrn",I4.4,".hdf5")')n
    inquire(file=trim(fname),exist=exist)
    if (exist) then
      call hdf5_read(fname,"/parameters","nwloc",nwloc)
      do iwloc=1,nwloc
        write(c1,'(I8.8)')iwloc
        call hdf5_read(fname,"/iwloc/"//c1,"iw",iw)
        call hdf5_read(fname,"/iwloc/"//c1,"w",w(iw))
        do it=1,ntr
          write(c2,'("t",I7.7)')it
          call hdf5_read(fname,"/iwloc/"//c1//"/"//c2,"u4",u4mtrx(1,1,it),&
            (/nwt,nwt/))
        enddo
        do i=1,nu4
          i1=iwtidx(iu4(1,i),iu4(2,i),iu4(3,i),iu4(4,i),iu4(5,i))
          i2=iwtidx(iu4(6,i),iu4(7,i),iu4(8,i),iu4(9,i),iu4(10,i))
          it=vtridx(iu4(11,i),iu4(12,i),iu4(13,i))
          u4(i,iw)=u4mtrx(i1,i2,it)*ha2ev
        enddo
      enddo !iwloc
    endif !exist
  enddo !n 
! write to file   
  open(150,file=trim(adjustl(fout)),status="REPLACE",form="FORMATTED")
!  write(150,'("# Wannier functions : ",100I4)')iwann_
!  write(150,'("#")')
!  write(150,'("# Screened U(w=0) matrix")')
!  write(150,'("# real part ")')  
!  do i=1,nwann_
!    write(150,'("# ",100F16.8)')(dreal(uscrn(i,j,1)),j=1,nwann_)
!  enddo
!  write(150,'("# image part ")')  
!  do i=1,nwann_
!    write(150,'("# ",100F16.8)')(dimag(uscrn(i,j,1)-zeps),j=1,nwann_)
!  enddo
!  write(150,'("#")')    
!  write(150,'("# Screened J(w=0) matrix")')
!  write(150,'("# real part ")')  
!  do i=1,nwann_
!    write(150,'("# ",100F16.8)')(dreal(jscrn(i,j,1)),j=1,nwann_)
!  enddo
!  write(150,'("# image part ")')  
!  do i=1,nwann_
!    write(150,'("# ",100F16.8)')(dimag(jscrn(i,j,1)-zeps),j=1,nwann_)
!  enddo
!  write(150,'("#")')    
!  write(150,'("# UJ matrix")')
!  write(150,'("# real part ")')    
!  do i=1,2*nwann_
!    write(150,'("# ",100F16.8)')(dreal(ujmtrx(i,j)),j=1,2*nwann_)
!  enddo
!  write(150,'("# image part ")')    
!  do i=1,2*nwann_
!    write(150,'("# ",100F16.8)')(dimag(ujmtrx(i,j)-zeps),j=1,2*nwann_)
!  enddo
!  write(150,'("#")')
!  write(150,'("# columns (units are eV) : ")')
!  write(150,'("#   1 : energy ")')
!  write(150,'("#   2 : Re (average U) ")')
!  write(150,'("#   3 : Im (average U) ")')
!  write(150,'("#   4 : Re (average J) ")')
!  write(150,'("#   5 : Im (average J) ")')
!  write(150,'("#   6 : Re (U_{11}) ")')
!  write(150,'("#   7 : Im (U_{11}) ")')
!  write(150,'("#   8 : Re (U_{22}) ")')
!  write(150,'("#   9 : Im (U_{22}) ")')
!  write(150,'("#  10 : ... ")')
!  write(150,'("#")')
  do iw=1,nw
    write(150,'(100F16.8)')w(iw)*ha2ev,(dreal(u4(i,iw)),&
      dimag(u4(i,iw)-zeps),i=1,nu4)
  enddo
  close(150)
endif !mode.eq.1

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
!          call hdf5_read(fname,"/iwloc/"//c8,"uscrn",uscrn(1),(/nwt/))
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
end subroutine

