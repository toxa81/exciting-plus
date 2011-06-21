subroutine gensmesh_cubed
use modmain
use mod_sic
implicit none
integer n,i,i1,i2,i3,ix,itp,lm,nf
real(8) px,py,pz,a,alpha,beta,t1,t2,t3,p1,p2,p3
logical texist
!

inquire(file="cubed.in",exist=texist)
if (texist) then
  open(177,file="cubed.in",status="OLD",form="FORMATTED")
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

!subroutine gen_bsht
!use modmain
!use mod_sic
!implicit none
!integer ierr
!real(8), allocatable :: o(:,:)
!real(8), allocatable :: a(:,:)
!complex(8), allocatable :: zo(:,:)
!complex(8), allocatable :: za(:,:)
!logical texist,tgen
!integer lmmaxwan_,s_ntp_,itp
!real(8), allocatable :: s_tp_(:,:)
!!
!tgen=.false.
!inquire(file="bsht",exist=texist)
!if (texist) then
!  if (mpi_grid_root()) then
!    open(300,file="bsht",form="unformatted")
!    read(300)lmmaxwan_
!    read(300)s_ntp_
!    allocate(s_tp_(2,s_ntp_))
!    read(300)s_tp_
!    if (lmmaxwan_.ne.lmmaxwan) tgen=.true.
!    if (s_ntp_.ne.s_ntp) then
!      tgen=.true.
!    else
!      do itp=1,s_ntp
!        if (s_tp_(1,itp).ne.s_tp(1,itp)) tgen=.true.
!        if (s_tp_(2,itp).ne.s_tp(2,itp)) tgen=.true.
!      enddo
!    endif
!    if (.not.tgen) then
!      read(300)s_ylmb
!      read(300)s_rlmb
!    endif
!    close(300)
!    deallocate(s_tp_)
!  endif
!  call mpi_grid_bcast(tgen)
!else
!  tgen=.true.
!endif
!if (tgen) then
!  if (mpi_grid_root()) write(*,'("[gen_bsht] generating transformation matrices")')
!  allocate(zo(lmmaxwan,lmmaxwan))
!  allocate(za(lmmaxwan,lmmaxwan))
!! compute overlap matrix o_{lm,l'm'}=<Y_{lm}|Y_{l'm'}> 
!!  note: zgemm computes conjugated overlap matrix
!  call zgemm('N','C',lmmaxwan,lmmaxwan,s_ntp,zone,s_ylmf,lmmaxwan,&
!    s_ylmf,lmmaxwan,zzero,zo,lmmaxwan)
!  zo=dconjg(zo)
!! calculate S=O^{-1/2}
!  call isqrtzhe(lmmaxwan,zo,ierr) 
!  if (ierr.ne.0) then
!    write(*,'("Error(gen_bsht): overlap matrix of complex spherical harmonics&
!     & is degenerate")')
!    call pstop
!  endif
!  call zgemm('N','C',lmmaxwan,lmmaxwan,lmmaxwan,zone,zo,lmmaxwan,zo,lmmaxwan,&
!    zzero,za,lmmaxwan)
!  call zgemm('C','T',s_ntp,lmmaxwan,lmmaxwan,zone,s_ylmf,lmmaxwan,za,lmmaxwan,&
!    zzero,s_ylmb,s_ntp)
!  deallocate(zo,za)
!  allocate(o(lmmaxwan,lmmaxwan))
!  allocate(a(lmmaxwan,lmmaxwan))
!! compute overlap matrix o_{lm,l'm'}=<R_{lm}|R_{l'm'}> 
!  call dgemm('N','T',lmmaxwan,lmmaxwan,s_ntp,1.d0,s_rlmf,lmmaxwan,&
!    s_rlmf,lmmaxwan,0.d0,o,lmmaxwan)
!! calculate S=O^{-1/2}
!  call isqrtdsy(lmmaxwan,o,ierr) 
!  if (ierr.ne.0) then
!    write(*,'("Error(gen_bsht): overlap matrix of real spherical harmonics&
!     & is degenerate")')
!    call pstop
!  endif
!  call dgemm('N','T',lmmaxwan,lmmaxwan,lmmaxwan,1.d0,o,lmmaxwan,o,lmmaxwan,&
!    0.d0,a,lmmaxwan)
!  call dgemm('T','T',s_ntp,lmmaxwan,lmmaxwan,1.d0,s_rlmf,lmmaxwan,a,lmmaxwan,&
!    0.d0,s_rlmb,s_ntp)
!   deallocate(a,o)
!else
!  call mpi_grid_bcast(s_ylmb(1,1),s_ntp*lmmaxwan)
!  call mpi_grid_bcast(s_rlmb(1,1),s_ntp*lmmaxwan)
!endif
!if (mpi_grid_root().and.tgen) then
!  open(300,file="bsht",form="unformatted",status="replace")
!  write(300)lmmaxwan
!  write(300)s_ntp
!  write(300)s_tp
!  write(300)s_ylmb
!  write(300)s_rlmb
!  close(300)
!endif
!return
!end


