subroutine genwu(ng,iw,chi0,vcgq,qnm)
use modmain
implicit none
integer, intent(in) :: ng
integer, intent(in) :: iw
complex(8), intent(in) :: chi0(ng,ng)
real(8), intent(in) :: vcgq(ng)
character(*), intent(in) :: qnm

complex(8), allocatable :: epsilon(:,:)
complex(8), allocatable :: mtrx1(:,:)
!complex(8), allocatable :: uscrn(:,:)
!complex(8), allocatable :: ubare(:,:)
integer ig1,ig2,n1,n2
integer it1,i,j
character*100 fname
integer, external :: hash
character*12 c12

allocate(epsilon(ng,ng))
allocate(mtrx1(ng,ng))

!! rpa kernel
!mtrx1=dcmplx(0.d0,0.d0)
!do i=1,ng
!  mtrx1(i,i)=vcgq(i)**2
!enddo
!! compute matrix epsilon=1-chi0*v
!epsilon=dcmplx(0.d0,0.d0)
!do i=1,ng
!  epsilon(i,i)=dcmplx(1.d0,0.d0)
!enddo
!call zgemm('N','N',ng,ng,ng,dcmplx(-1.d0,0.d0),chi0,ng,mtrx1,ng,&
!  dcmplx(1.d0,0.d0),epsilon,ng)
!! invert epsilon matrix
!call invzge(epsilon,ng)
!! compute chi=epsilon^-1 * chi0
!call zgemm('N','N',ng,ng,ng,dcmplx(1.d0,0.d0),epsilon,ng,chi0,ng,&
!  dcmplx(0.d0,0.d0),mtrx1,ng)

! compute screened Coulomb potential using "symmetrized" dielectric function
do ig1=1,ng
  do ig2=1,ng
    epsilon(ig1,ig2)=-vcgq(ig1)*chi0(ig1,ig2)*vcgq(ig2)
  enddo
  epsilon(ig1,ig1)=zone+epsilon(ig1,ig1)
enddo
call invzge(epsilon,ng)
do ig1=1,ng
  do ig2=1,ng
    mtrx1(ig1,ig2)=vcgq(ig1)*epsilon(ig1,ig2)*vcgq(ig2)
  enddo
enddo

do i=1,ntrmegqwan
  if (itrmegqwan(1,i).eq.0.and.itrmegqwan(2,i).eq.0.and.&
      itrmegqwan(3,i).eq.0) then
    it1=i
    exit
  endif
enddo
!allocate(uscrn(nwann,nwann))
!allocate(ubare(nwann,nwann))
!uscrn=zzero
! compute screened u
do ig1=1,ng
  do ig2=1,ng
    do n1=1,nwann
      do n2=1,nwann
        uscrnwan(n1,n2,iw)=uscrnwan(n1,n2,iw)+dconjg(megqwan(imegqwan(n1,n1),it1,ig1))*&
          mtrx1(ig1,ig2)*megqwan(imegqwan(n2,n2),it1,ig2)
      enddo
    enddo
  enddo
enddo
! compute bare u
!ubare=zzero
if (iw.eq.1) then
  do ig1=1,ng
    do n1=1,nwann
      do n2=1,nwann
        ubarewan(n1,n2)=ubarewan(n1,n2)+dconjg(megqwan(imegqwan(n1,n1),it1,ig1))*&
          (vcgq(ig1)**2)*megqwan(imegqwan(n2,n2),it1,ig1)
      enddo
    enddo
  enddo
endif

! write block of W matrix
!if (ng.gt.10) then
!  n1=10
!else
!  n1=ng
!endif
!fname=trim(qnm)//"_W__.txt"
!open(170,file=trim(fname),status='replace',form='formatted')
!write(170,'("Screened W matrix")')
!write(170,'("real part")')
!do ig1=1,n1
!  write(170,'(100F12.6)')(dreal(mtrx1(ig1,ig2)),ig2=1,n1)
!enddo
!write(170,'("imag part")')
!do ig1=1,n1
!  write(170,'(100F12.6)')(dimag(mtrx1(ig1,ig2)),ig2=1,n1)
!enddo
!close(170)
!
!fname=trim(qnm)//"_U__.txt"
!open(170,file=trim(fname),status='replace',form='formatted')
!write(170,'("Screened U matrix")')
!write(170,'("real part")')
!do i=1,nwann
!  write(170,'(100F12.6)')(dreal(uscrn(i,j)),j=1,nwann)
!enddo
!write(170,'("imag part")')
!do i=1,nwann
!  write(170,'(100F12.6)')(dimag(uscrn(i,j)),j=1,nwann)
!enddo
!write(170,*)
!write(170,'("Bare U matrix")')
!write(170,'("real part")')
!do i=1,nwann
!  write(170,'(100F12.6)')(dreal(ubare(i,j)),j=1,nwann)
!enddo
!write(170,'("imag part")')
!do i=1,nwann
!  write(170,'(100F12.6)')(dimag(ubare(i,j)),j=1,nwann)
!enddo
!
!fname=trim(qnm)//"_U"
!open(170,file=trim(fname),status='replace',form='unformatted')
!write(170)uscrn,ubare
!close(170)

!fname=trim(qnm)//"_lr.hdf5"
!if (mpi_grid_root(dims=(/dim_k,dim_b/))) then
!  write(c12,'("/iw/",I8.8)')iw
!  call write_real8(lr_w(iw),2,trim(fname),c12,'w')
!  call write_real8_array(uscrn,3,(/2,nwann,nwann/),trim(fname),c12,'uscrn')
!endif

deallocate(epsilon,mtrx1)!,uscrn,ubare)

return 
end