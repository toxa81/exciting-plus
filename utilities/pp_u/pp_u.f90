program pp_u
use mod_hdf5
implicit none

character*100 fname
integer nw,size,ntr_uscrn,it,n,nwann,iw,iwloc,nwloc,i,j
integer n1,n2
real(8), allocatable :: w(:)
complex(8), allocatable :: uscrn(:)
complex(8), allocatable :: jscrn(:)

complex(8), allocatable :: uscrn_(:,:,:)
complex(8), allocatable :: jscrn_(:,:,:)
complex(8), allocatable :: uscrn2_(:,:)
complex(8), allocatable :: ujmtrx(:,:)
complex(8) zf(2)
real(8) t1,uavg
character*10 c1,c2,c3,c4
character*8 c8
logical exist
real(8), parameter :: ha2ev = 27.21138386d0
integer vtl_(3)
integer nwann_,is1,is2
integer, allocatable :: iwann_(:)
integer nlist_
integer, allocatable :: ilist_(:)
integer nmegqwan
integer, allocatable :: imegqwan(:,:)
integer, allocatable :: iwm(:,:)
character*100 str,fout
integer mode
integer ngq,ngridk(3)

call hdf5_initialize
open(150,file="pp_u.in",form="FORMATTED",status="OLD")
read(150,'(A)')str
if (trim(adjustl(str)).eq."onsite") then
  mode=0
  vtl_=0
  read(150,*)nwann_
  allocate(iwann_(nwann_))
  read(150,*)iwann_
  close(150)
else if (trim(adjustl(str)).eq."offsite") then
  mode=1
else if (trim(adjustl(str)).eq."list") then
  mode=2
  read(150,*)nlist_
  allocate(ilist_(nlist_))
  read(150,*)ilist_
endif
write(*,*)trim(adjustl(str))
write(*,*)mode

fname="uscrn0000.hdf"
inquire(file=trim(fname),exist=exist)
if (exist) then
  call hdf5_read(fname,"/parameters","nw",nw)
  call hdf5_read(fname,"/parameters","nwann",nwann)
  call hdf5_read(fname,"/parameters","size",size)
  call hdf5_read(fname,"/parameters","nmegqwan",nmegqwan)
  allocate(imegqwan(5,nmegqwan))
  call hdf5_read(fname,"/parameters","imegqwan",imegqwan(1,1),(/5,nmegqwan/))  
  call hdf5_read(fname,"/parameters","ngq",ngq)
  call hdf5_read(fname,"/parameters","ngridk",ngridk(1),(/3/))  
  allocate(w(nw))
  allocate(uscrn(nmegqwan))
  allocate(jscrn(nmegqwan))
  write(c1,'(I6)')ngridk(1)
  write(c2,'(I6)')ngridk(2)
  write(c3,'(I6)')ngridk(3)
  write(c4,'(I8)')ngq
  fout="cRPA__"//trim(adjustl(c1))//"x"//trim(adjustl(c2))//"x"//&
    trim(adjustl(c3))//"_kgrid__"//trim(adjustl(c4))//"_Gv__.dat"
  if (mode.eq.0) then
    allocate(iwm(nwann_,nwann_)) 
    allocate(uscrn_(nwann_,nwann_,nw))
    allocate(jscrn_(nwann_,nwann_,nw))
    uscrn_=dcmplx(0.d0,0.d0)
    jscrn_=dcmplx(0.d0,0.d0)
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
          call hdf5_read(fname,"/iwloc/"//c8,"jscrn",jscrn(1),(/nmegqwan/))
          do n1=1,nwann_
            do n2=1,nwann_
              uscrn_(n1,n2,iw)=uscrn(iwm(n1,n2))*ha2ev
              jscrn_(n1,n2,iw)=jscrn(iwm(n1,n2))*ha2ev
            enddo
          enddo
        enddo
      endif
    enddo
    
    allocate(ujmtrx(2*nwann_,2*nwann_))
    do is1=1,2
      do is2=1,2
        do i=1,nwann_
          do j=1,nwann_
	    if (is1.eq.is2) then
	      ujmtrx((is1-1)*nwann_+i,(is2-1)*nwann_+j)=uscrn_(i,j,1)-jscrn_(i,j,1)
	    else
	      ujmtrx((is1-1)*nwann_+i,(is2-1)*nwann_+j)=uscrn_(i,j,1)
	    endif
	  enddo
	enddo
      enddo
    enddo
    do i=1,2*nwann_
      write(*,'(255F12.6)')(dreal(ujmtrx(i,j)),j=1,2*nwann_)
    enddo
    deallocate(ujmtrx)
    
    open(150,file=trim(adjustl(fout)),status="REPLACE",form="FORMATTED")
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
    write(150,'("# Screened J(w=0) matrix")')
    write(150,'("#  real part")')
    do i=1,nwann_
      write(150,'("# ",100F12.6)')(dreal(jscrn_(i,j,1)),j=1,nwann_)
    enddo
    write(150,'("#  imag part")')
    do i=1,nwann_
      write(150,'("# ",100F12.6)')(dimag(jscrn_(i,j,1)),j=1,nwann_)
    enddo
    write(150,'("#")')
    write(150,'("# columns : ")')
    write(150,'("#   1 : energy ")')
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
  endif !mode.eq.0
  
  if (mode.eq.2) then
    allocate(uscrn2_(nlist_,nw))
    uscrn2_=dcmplx(0.d0,0.d0)
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
          do i=1,nlist_
            uscrn2_(i,iw)=uscrn(ilist_(i))*ha2ev
          enddo
        enddo
      endif
    enddo
    
    open(150,file=trim(adjustl(fout)),status="REPLACE",form="FORMATTED")
    write(150,'("#")')
    write(150,'("# columns : ")')
    write(150,'("#   1 : energy ")')
    write(150,'("#   2 : Re(U(first transition in list)) ")')
    write(150,'("#   3 : Im(U(first transition in list)) ")')
    write(150,'("#   4 : Re(U(second transition in list)) ")')
    write(150,'("#   5 : Im(U(second transition in list)) ")')    
    write(150,'("#   5 : ... ")')
    write(150,'("#")')
    do iw=1,nw
      write(150,'(100G18.10)')w(iw)*ha2ev,(dreal(uscrn2_(i,iw)),&
        dimag(uscrn2_(i,iw)),i=1,nlist_)
    enddo
    close(150)
  endif
endif
call hdf5_finalize
end

subroutine pstop
stop
end