subroutine gensmesh_cubed
use modmain
use mod_sic
implicit none
integer n,i,i1,i2,i3,nx,ix,itp,lm,nf
real(8) px,py,pz,a,alpha,beta,t1,t2,t3,p1,p2,p3
real(8) x0(3),x1(3),x2(3),x3(3)
real(8), allocatable :: x(:,:)
integer, allocatable :: facets(:,:)

n=10

nx=2*n*n+2*n*(n-2)+2*(n-2)*(n-2)
allocate(x(3,nx))

ix=0
do i1=1,n
  do i2=1,n
    alpha=-pi/4.d0+dble(i1-1)/(n-1)*pi/2
    beta=-pi/4.d0+dble(i2-1)/(n-1)*pi/2
    px=tan(alpha)
    py=tan(beta)
    ix=ix+2
    x(:,ix-1)=(/px,py,-1.d0/)
    x(:,ix)=(/px,py,1.d0/)
  enddo
enddo
do i1=1,n
  do i2=2,n-1
    alpha=-pi/4.d0+dble(i1-1)/(n-1)*pi/2
    beta=-pi/4.d0+dble(i2-1)/(n-1)*pi/2
    px=tan(alpha)
    pz=tan(beta)
    ix=ix+2
    x(:,ix-1)=(/px,-1.d0,pz/)
    x(:,ix)=(/px,1.d0,pz/)
  enddo
enddo
do i1=2,n-1
  do i2=2,n-1
    alpha=-pi/4.d0+dble(i1-1)/(n-1)*pi/2
    beta=-pi/4.d0+dble(i2-1)/(n-1)*pi/2
    py=tan(alpha)
    pz=tan(beta)
    ix=ix+2
    x(:,ix-1)=(/-1.d0,py,pz/)
    x(:,ix)=(/1.d0,py,pz/)
  enddo
enddo
do i=1,nx
  x(:,i)=x(:,i)/sqrt(sum(x(:,i)**2))
enddo

!allocate(facets(3,nx*5))
!call find_facets(nx,x,nf,facets)

!s_ntp=nx
!if (allocated(s_tpw)) deallocate(s_tpw)
!allocate(s_tpw(s_ntp))
!s_tpw=0.d0

!do i=1,nf
!  x1=x(:,facets(1,i)) 
!  x2=x(:,facets(2,i)) 
!  x3=x(:,facets(3,i))
!  call r3cross(x2,x3,x0)
!  t2=dot_product(x1,x0)
!  p1=sqrt(sum(x1**2))
!  p2=sqrt(sum(x2**2))
!  p3=sqrt(sum(x3**2))
!  t3=p1*p2*p3+dot_product(x1,x2)*p3+dot_product(x1,x3)*p2+&
!    dot_product(x2,x3)*p1
!  t1=atan(t2/t3)*2.d0
!  !write(*,*)t1
!  s_tpw(facets(1,i))=s_tpw(facets(1,i))+t1/3.d0
!  s_tpw(facets(2,i))=s_tpw(facets(2,i))+t1/3.d0
!  s_tpw(facets(3,i))=s_tpw(facets(3,i))+t1/3.d0
!enddo
!deallocate(facets)

s_ntp=nx
if (allocated(s_tp)) deallocate(s_tp)
allocate(s_tp(2,s_ntp))
if (allocated(s_x)) deallocate(s_x)
allocate(s_x(3,s_ntp))
if (allocated(s_tpw)) deallocate(s_tpw)
allocate(s_tpw(s_ntp))
if (allocated(s_rlmf)) deallocate(s_rlmf)
allocate(s_rlmf(lmmaxwan,s_ntp))
if (allocated(s_ylmf)) deallocate(s_ylmf)
allocate(s_ylmf(lmmaxwan,s_ntp))
if (allocated(s_rlmb)) deallocate(s_rlmb)
allocate(s_rlmb(s_ntp,lmmaxwan))
if (allocated(s_ylmb)) deallocate(s_ylmb)
allocate(s_ylmb(s_ntp,lmmaxwan))
s_x(:,1:s_ntp)=x(:,1:nx)
do itp=1,s_ntp
  call sphcrd(s_x(1,itp),a,s_tp(1,itp))
  s_tpw(itp)=fourpi/s_ntp
enddo
do itp=1,s_ntp
  call genrlm(lmaxwan,s_tp(1,itp),s_rlmf(1,itp))
  call genylm(lmaxwan,s_tp(1,itp),s_ylmf(1,itp))
  do lm=1,lmmaxwan
    s_rlmb(itp,lm)=s_rlmf(lm,itp)*s_tpw(itp)
    s_ylmb(itp,lm)=dconjg(s_ylmf(lm,itp))*s_tpw(itp) 
  enddo
enddo
call gen_bsht
!do itp=1,2
!  do lm=1,5
!    write(*,*)"lm=",lm,s_rlmb(itp,lm)/s_rlmf(lm,itp)
!  enddo
!enddo
!call bstop
!open(200,file="sph_x.dat",form="formatted",status="replace")
!do i=1,s_ntp
!  write(200,'(3(G18.10,","),G18.10)')s_x(:,i),1.d0
!enddo
!close(200)
!open(200,file="sph_tp.dat",form="formatted",status="replace")
!do i=1,s_ntp
!  write(200,'(G18.10,",",G18.10)')s_tp(:,i)
!enddo
!close(200)
return
end subroutine


subroutine gen_bsht
use modmain
use mod_sic
implicit none
integer ierr
real(8), allocatable :: o(:,:)
real(8), allocatable :: a(:,:)
complex(8), allocatable :: zo(:,:)
complex(8), allocatable :: za(:,:)
logical texist,tgen
integer lmmaxwan_,s_ntp_,itp
real(8), allocatable :: s_tp_(:,:)
!
tgen=.false.
inquire(file="bsht",exist=texist)
if (texist) then
  if (mpi_grid_root()) then
    open(300,file="bsht",form="unformatted")
    read(300)lmmaxwan_
    read(300)s_ntp_
    allocate(s_tp_(2,s_ntp_))
    read(300)s_tp_
    if (lmmaxwan_.ne.lmmaxwan) tgen=.true.
    if (s_ntp_.ne.s_ntp) then
      tgen=.true.
    else
      do itp=1,s_ntp
        if (s_tp_(1,itp).ne.s_tp(1,itp)) tgen=.true.
        if (s_tp_(2,itp).ne.s_tp(2,itp)) tgen=.true.
      enddo
    endif
    if (.not.tgen) then
      read(300)s_ylmb
      read(300)s_rlmb
    endif
    close(300)
    deallocate(s_tp_)
  endif
  call mpi_grid_bcast(tgen)
else
  tgen=.true.
endif
if (tgen) then
  if (mpi_grid_root()) write(*,'("[gen_bsht] generating transformation matrices")')
  allocate(zo(lmmaxwan,lmmaxwan))
  allocate(za(lmmaxwan,lmmaxwan))
! compute overlap matrix o_{lm,l'm'}=<Y_{lm}|Y_{l'm'}> 
!  note: zgemm computes conjugated overlap matrix
  call zgemm('N','C',lmmaxwan,lmmaxwan,s_ntp,zone,s_ylmf,lmmaxwan,&
    s_ylmf,lmmaxwan,zzero,zo,lmmaxwan)
  zo=dconjg(zo)
! calculate S=O^{-1/2}
  call isqrtzhe(lmmaxwan,zo,ierr) 
  if (ierr.ne.0) then
    write(*,'("Error(gen_bsht): overlap matrix of complex spherical harmonics&
     & is degenerate")')
    call pstop
  endif
  call zgemm('N','C',lmmaxwan,lmmaxwan,lmmaxwan,zone,zo,lmmaxwan,zo,lmmaxwan,&
    zzero,za,lmmaxwan)
  call zgemm('C','T',s_ntp,lmmaxwan,lmmaxwan,zone,s_ylmf,lmmaxwan,za,lmmaxwan,&
    zzero,s_ylmb,s_ntp)
  deallocate(zo,za)
  allocate(o(lmmaxwan,lmmaxwan))
  allocate(a(lmmaxwan,lmmaxwan))
! compute overlap matrix o_{lm,l'm'}=<R_{lm}|R_{l'm'}> 
  call dgemm('N','T',lmmaxwan,lmmaxwan,s_ntp,1.d0,s_rlmf,lmmaxwan,&
    s_rlmf,lmmaxwan,0.d0,o,lmmaxwan)
! calculate S=O^{-1/2}
  call isqrtdsy(lmmaxwan,o,ierr) 
  if (ierr.ne.0) then
    write(*,'("Error(gen_bsht): overlap matrix of real spherical harmonics&
     & is degenerate")')
    call pstop
  endif
  call dgemm('N','T',lmmaxwan,lmmaxwan,lmmaxwan,1.d0,o,lmmaxwan,o,lmmaxwan,&
    0.d0,a,lmmaxwan)
  call dgemm('T','T',s_ntp,lmmaxwan,lmmaxwan,1.d0,s_rlmf,lmmaxwan,a,lmmaxwan,&
    0.d0,s_rlmb,s_ntp)
   deallocate(a,o)
else
  call mpi_grid_bcast(s_ylmb(1,1),s_ntp*lmmaxwan)
  call mpi_grid_bcast(s_rlmb(1,1),s_ntp*lmmaxwan)
endif
if (mpi_grid_root().and.tgen) then
  open(300,file="bsht",form="unformatted",status="replace")
  write(300)lmmaxwan
  write(300)s_ntp
  write(300)s_tp
  write(300)s_ylmb
  write(300)s_rlmb
  close(300)
endif
return
end


