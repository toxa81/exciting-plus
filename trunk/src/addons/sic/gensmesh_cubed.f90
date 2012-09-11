subroutine gensmesh_cubed
use modmain
use mod_sic
implicit none
integer n,i1,i2,ix,itp
real(8) px,py,pz,a,alpha,beta
logical texist
!
inquire(file="cubed.in",exist=texist)
if (texist) then
  open(177,file="cubed.in",status="OLD",form="FORMATTED")
  read(177,*)lmaxwan
  lmmaxwan=(lmaxwan+1)**2
  read(177,*)s_ntp
  if (mpi_grid_root()) then
    write(*,'("[gensmesh_cubed] reading points from file")')
    write(*,'("[gensmesh_cubed] number of points : ",I4)')s_ntp
  endif
  if (allocated(s_tp)) deallocate(s_tp)
  allocate(s_tp(2,s_ntp))
  if (allocated(s_x)) deallocate(s_x)
  allocate(s_x(3,s_ntp))
  if (allocated(s_tpw)) deallocate(s_tpw)
  allocate(s_tpw(s_ntp))
  do itp=1,s_ntp
    read(177,*)s_x(1,itp),s_x(2,itp),s_x(3,itp),s_tpw(itp)
  enddo
  close(177)
  do itp=1,s_ntp
    call sphcrd(s_x(1,itp),a,s_tp(1,itp))
  enddo
else
  n=sic_smesh_n
  s_ntp=2*n*n+2*n*(n-2)+2*(n-2)*(n-2)
  if (mpi_grid_root()) then
    write(*,'("[gensmesh_cubed] number of points : ",I4)')s_ntp
  endif
  if (allocated(s_tp)) deallocate(s_tp)
  allocate(s_tp(2,s_ntp))
  if (allocated(s_x)) deallocate(s_x)
  allocate(s_x(3,s_ntp))
  if (allocated(s_tpw)) deallocate(s_tpw)
  allocate(s_tpw(s_ntp))
  ix=0
  do i1=1,n
    do i2=1,n
      alpha=-pi/4.d0+dble(i1-1)/(n-1)*pi/2
      beta=-pi/4.d0+dble(i2-1)/(n-1)*pi/2
      px=tan(alpha)
      py=tan(beta)
      ix=ix+2
      s_x(:,ix-1)=(/px,py,-1.d0/)
      s_x(:,ix)=(/px,py,1.d0/)
    enddo
  enddo
  do i1=1,n
    do i2=2,n-1
      alpha=-pi/4.d0+dble(i1-1)/(n-1)*pi/2
      beta=-pi/4.d0+dble(i2-1)/(n-1)*pi/2
      px=tan(alpha)
      pz=tan(beta)
      ix=ix+2
      s_x(:,ix-1)=(/px,-1.d0,pz/)
      s_x(:,ix)=(/px,1.d0,pz/)
    enddo
  enddo
  do i1=2,n-1
    do i2=2,n-1
      alpha=-pi/4.d0+dble(i1-1)/(n-1)*pi/2
      beta=-pi/4.d0+dble(i2-1)/(n-1)*pi/2
      py=tan(alpha)
      pz=tan(beta)
      ix=ix+2
      s_x(:,ix-1)=(/-1.d0,py,pz/)
      s_x(:,ix)=(/1.d0,py,pz/)
    enddo
  enddo
  do itp=1,s_ntp
    call sphcrd(s_x(1,itp),a,s_tp(1,itp))
    s_x(:,itp)=s_x(:,itp)/a
    s_tpw(itp)=fourpi/s_ntp
  enddo
endif
return
end subroutine
