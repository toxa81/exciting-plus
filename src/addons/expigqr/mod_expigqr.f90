module mod_expigqr
use mod_wannier
implicit none

! dimension for interband transitions
!integer, parameter :: dim_b=3

! number of G-vectors for matrix elements
integer ngvecme

! if wave-function is a 2-component spinor then e^{-i(G+q)x} must be 
! a 2x2 matrix in spin space; expigqr22 controls the valuse of this 2x2 matrix
! expigqr22=1 : diagonal matrix for charge response
integer expigqr22
data expigqr22/1/

! total number of matrix elements <nk|e^{-i(G+q)x}|n'k+q> in the Bloch basis
! for a given local k-point
integer, allocatable :: nmegqblh(:)

! maximum total number of matrix elements <nk|e^{-i(G+q)x}|n'k+q> 
! over all k-points
!integer nmegqblhtotmax

! local number of interband transitions and offset (each processor along 
! dim_b does it's local fraction of nmegqblh(ikloc) transitions)
!   1-st index: 1: local number of transitions
!               2: offset to compute global index
!integer, allocatable :: nmegqblh(:,:)

! maximum local number of matrix elements <nk|e^{-i(G+q)x}|n'k+q> 
! over all k-points
!integer nmegqblhlocmax

! bands (n,n') for matrix elements <nk|e^{-i(G+q)x}|n'k+q>  
!   1-st index :  1: n at k
!                 2: n' at k+q
!   2-nd index : global index of pair of bands (n,n')
!   3-rd index : k-point
integer, allocatable :: bmegqblh(:,:,:)

! matrix elements <nk|e^{-i(G+q)x}|n'k+q> in the Bloch basis
!   1-st index : local index of pair of bands (n,n')
!   2-nd index : G-vector
!   3-rd index : k-point
complex(8), allocatable :: megqblh(:,:,:)

! adjoint matrix elements <n,k-q|e^{-i(G+q)x}|n'k> in the Bloch basis
complex(8), allocatable :: amegqblh(:,:,:)
! number of adjoint matrix elements 
integer, allocatable :: namegqblh(:)
! band indices of adjoint matrix elements
integer, allocatable :: bamegqblh(:,:,:)

! interval of bands to take for matrix elements <nk|e^{-i(G+q)x}|n'k+q>
real(8) megq_include_bands(2)
data megq_include_bands/-100.1d0,100.1d0/

! minimum interband transition energy
real(8) lr_min_e12

! low level switch: compute matrix elements of e^{i(G+q)x} in the basis of
!   Wannier functions; depends on crpa and wannier_chi0_chi
logical wannier_megq
data wannier_megq/.false./

! minimum and maximum cutoff values for matrix elements in Wannier basis
real(8) megqwan_cutoff(2)
data megqwan_cutoff/-0.0d0,1000.d0/

real(8) megqwan_mindist
data megqwan_mindist/-0.0d0/
real(8) megqwan_maxdist
data megqwan_maxdist/0.1d0/

complex(8), allocatable :: megqwan(:,:)

integer nwann_include
data nwann_include/0/
integer, allocatable :: iwann_include(:)

integer nmegqblhwanmax
integer, allocatable :: nmegqblhwan(:)
integer, allocatable :: imegqblhwan(:,:)

complex(8), allocatable :: wann_c_jk(:,:,:)

integer ngntujumax
integer, allocatable :: ngntuju(:,:)
integer(2), allocatable :: igntuju(:,:,:,:)
complex(8), allocatable :: gntuju(:,:,:)

! array for k+q points
!  1-st index: index of k-point in BZ
!  2-nd index: 1: index of k'=k+q-K
!              2: index of K-vector which brings k+q to first BZ
!              3: index of k'=k-q point
integer, allocatable :: idxkq(:,:)

type(wannier_transitions) :: megqwantran

end module
