subroutine disentangle(evalsv_,wann_c_,evecsv_)
use modmain
use mod_wannier
implicit none
real(8), intent(inout) :: evalsv_(nstsv)
complex(8), intent(inout) :: wann_c_(nwantot,nstsv)
complex(8), intent(inout) :: evecsv_(nstsv,nstsv)
complex(8), allocatable :: tmtrx1(:,:)
complex(8), allocatable :: tmtrx2(:,:)
complex(8), allocatable :: zv1(:)
complex(8), allocatable :: hwann(:,:)
complex(8), allocatable :: hrest(:,:)
real(8), allocatable :: ewann(:)
real(8), allocatable :: erest(:)
complex(8), allocatable :: wann_c_new(:,:)
complex(8), allocatable :: evecsv_new(:,:)
integer n,n1,m1,m2,i,j
complex(8) zt1
real(8) t1

! check initial orthogonality of eigen-vectors
t1=0.d0
do m1=1,nstsv
  do m2=1,nstsv
    zt1=zzero
    do j=1,nstsv
      zt1=zt1+dconjg(evecsv_(j,m1))*evecsv_(j,m2)
    enddo
    if (m1.eq.m2) zt1=zt1-zone
    t1=max(abs(zt1),t1)
  enddo
enddo
if (t1.gt.1d-12) then
  write(*,'("Error(disentangle): original eigen-vectors are not orthonormal")')
  write(*,*)"  maximum deviation : ",t1
  call pstop
endif

! check orthogonality of Wannier states
t1=0.d0
do m1=1,nwantot
  do m2=1,nwantot
    zt1=zzero
    do j=1,nstsv
      zt1=zt1+dconjg(wann_c_(m1,j))*wann_c_(m2,j)
    enddo
    if (m1.eq.m2) zt1=zt1-zone
    t1=max(abs(zt1),t1)
  enddo
enddo
if (t1.gt.1d-8) then
  write(*,'("Error(disentangle): Wannier states are not orthonormal")')
  write(*,*)"  maximum deviation : ",t1
  call pstop
endif


! transfomation matrix from eigen-vectors of original Hamiltonian to new basis
!  \tilde u_n=\sum_j tmtrx1(j,n) \psi_j
allocate(tmtrx1(nstsv,nstsv))
tmtrx1=zzero
! we already know first nwantot functions of a new basis - this are Bloch-sums
!   of Wannier functions
do n=1,nwantot
  tmtrx1(:,n)=wann_c_(n,:)
enddo
allocate(zv1(nstsv))
n1=nwantot+1
m1=1
do while (n1.le.nstsv)
  if (m1.gt.nstsv) then
    write(*,'("Error(disentangle): Gram-Schmidt orthonormalization failed")')
    call pstop
  endif
! trial vector
  tmtrx1(m1,n1)=zone
  zv1(:)=tmtrx1(:,n1)
! use Gram-Schmidt algorithm to orthogonalize trial vector
  do i=1,n1-1
    zv1(:)=zv1(:)-dot_product(tmtrx1(:,i),tmtrx1(:,n1))*tmtrx1(:,i)
  enddo
! vector norm
  zt1=dot_product(zv1,zv1)
  if (abs(zt1).gt.0.1d0) then
! normalize and add to the basis
    tmtrx1(:,n1)=zv1(:)/sqrt(abs(zt1))
    n1=n1+1
  endif
  m1=m1+1
enddo

! check orthogonality
t1=0.d0
do m1=1,nstsv
  do m2=1,nstsv
    zt1=zzero
    do j=1,nstsv
      zt1=zt1+dconjg(tmtrx1(j,m1))*tmtrx1(j,m2)
    enddo
    if (m1.eq.m2) zt1=zt1-zone
    t1=max(abs(zt1),t1)
  enddo
enddo
if (t1.gt.1d-8) then
  write(*,'("Error(disentangle): new basis is not orthonormal")')
  write(*,*)"  maximum deviation : ",t1
  call pstop
endif

allocate(hwann(nwantot,nwantot))
allocate(ewann(nwantot))
allocate(hrest(nstsv-nwantot,nstsv-nwantot))
allocate(erest(nstsv-nwantot))

hwann=zzero
do m1=1,nwantot
  do m2=1,nwantot
    do j=1,nstsv
      hwann(m1,m2)=hwann(m1,m2)+dconjg(tmtrx1(j,m1))*tmtrx1(j,m2)*evalsv_(j)
    enddo
  enddo
enddo

hrest=zzero
do m1=1,nstsv-nwantot
  do m2=1,nstsv-nwantot
    do j=1,nstsv
      hrest(m1,m2)=hrest(m1,m2)+dconjg(tmtrx1(j,m1+nwantot))*tmtrx1(j,m2+nwantot)*evalsv_(j)
    enddo
  enddo
enddo

call diagzhe(nwantot,hwann,ewann)
call diagzhe(nstsv-nwantot,hrest,erest)
evalsv_(1:nwantot)=ewann(1:nwantot)
evalsv_(nwantot+1:nstsv)=erest(1:nstsv-nwantot)

allocate(tmtrx2(nstsv,nstsv))
tmtrx2=zzero
! transformation matrix from psi_i of original Hamiltonian to \tilde \psi_i
!   of new Hamiltonian
! \tilde psi_i = \sum_j tmtrx2(j,i) \psi_j
do j=1,nstsv
  do i=1,nwantot
    do n=1,nwantot
      tmtrx2(j,i)=tmtrx2(j,i)+hwann(n,i)*tmtrx1(j,n)
    enddo
  enddo
  do i=1,nstsv-nwantot
    do n=1,nstsv-nwantot
      tmtrx2(j,i+nwantot)=tmtrx2(j,i+nwantot)+hrest(n,i)*tmtrx1(j,n+nwantot)
    enddo
  enddo
enddo
allocate(evecsv_new(nstsv,nstsv))
evecsv_new=zzero
do n=1,nstsv
  do i=1,nstsv
    do j=1,nstsv
      evecsv_new(n,i)=evecsv_new(n,i)+tmtrx2(j,i)*evecsv_(n,j)
    enddo
  enddo
enddo
evecsv_=evecsv_new

call invzge(tmtrx2,nstsv)

allocate(wann_c_new(nwantot,nstsv))
wann_c_new=zzero

do n=1,nwantot
  do i=1,nstsv
    do j=1,nstsv
      wann_c_new(n,i)=wann_c_new(n,i)+wann_c_(n,j)*tmtrx2(i,j)
    enddo
  enddo
enddo
wann_c_=wann_c_new

deallocate(tmtrx1)
deallocate(zv1)
deallocate(hwann)
deallocate(ewann)
deallocate(hrest)
deallocate(erest)
deallocate(tmtrx2)
deallocate(evecsv_new)

return
end
