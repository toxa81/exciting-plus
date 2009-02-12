program chi0diff
implicit none

integer gvecme1
integer gvecme2
integer ngvecchi
! number of G-vectors for matrix elements
integer ngvecme
! number of energy-mesh points
integer nepts
! energy mesh
complex(8), allocatable :: w(:)
! q-vector in lattice coordinates
real(8) vq0l(3)
! q-vector in Cartesian coordinates
real(8) vq0c(3)
! reduced q-vector in lattice coordinates
real(8) vq0rl(3)
! reduced q-vector in Cartesian coordinates
real(8) vq0rc(3)
! index of G-vector which brings q to first BZ
integer igq0
integer igq0s
! Kohn-Sham polarisability
complex(8), allocatable :: chi1(:,:,:,:)
complex(8), allocatable :: chi2(:,:,:,:)
! true polarisability
complex(8), allocatable :: chi(:,:,:)
! G+q vectors in Cartesian coordinates
real(8), allocatable :: vgq0c(:,:)
! length of G+q vectors
real(8), allocatable :: gq0(:)
! Coulomb potential 
real(8), allocatable :: vc(:)
! Fourier-transform of Ixc(r)=Bxc(r)/m(r)
complex(8), allocatable :: ixcft(:)

complex(8), allocatable :: epsilon_GqGq(:)
complex(8), allocatable :: chi_scalar(:)
complex(8), allocatable :: epsilon_eff(:)
complex(8), allocatable :: mtrx(:,:)
integer nspin_chi0
integer nspin_chi
integer gshme1,gshme2,spin_me

! allocatable arrays
real(8), allocatable :: func(:,:)
complex(8), allocatable :: epsilon(:,:)
integer, allocatable :: ipiv(:)

integer ie,ig,ngsh_me_,info,i,j,ig1,ig2,ispn
character*100 fname
integer iv(3)
real(8) diff
  
open(160,file='CHI1.OUT',form='unformatted',status='old')
read(160)nepts,igq0
read(160)gshme1,gshme2,gvecme1,gvecme2,ngvecme
allocate(w(nepts))
read(160)w(1:nepts)
read(160)vq0l(1:3)
read(160)vq0rl(1:3)
read(160)vq0c(1:3)
read(160)vq0rc(1:3)
read(160)spin_me,nspin_chi0
allocate(chi1(ngvecme,ngvecme,nepts,nspin_chi0))
do ie=1,nepts
  read(160)chi1(1:ngvecme,1:ngvecme,ie,1:nspin_chi0)
enddo
close(160)

open(160,file='CHI2.OUT',form='unformatted',status='old')
read(160)nepts,igq0
read(160)gshme1,gshme2,gvecme1,gvecme2,ngvecme
read(160)w(1:nepts)
read(160)vq0l(1:3)
read(160)vq0rl(1:3)
read(160)vq0c(1:3)
read(160)vq0rc(1:3)
read(160)spin_me,nspin_chi0
allocate(chi2(ngvecme,ngvecme,nepts,nspin_chi0))
do ie=1,nepts
  read(160)chi2(1:ngvecme,1:ngvecme,ie,1:nspin_chi0)
enddo
close(160)

do ie=1,nepts
  diff=sum(abs(chi1(:,:,ie,:)-chi2(:,:,ie,:)))
  write(*,*)'ie=',ie,'diff=',diff
enddo
diff=sum(abs(chi1(:,:,:,:)-chi2(:,:,:,:)))
write(*,*)'total diff=',diff
end
