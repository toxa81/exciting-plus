module modqpt
!-------------------------------!
!     q-point set variables     !
!-------------------------------!
! q-point grid sizes
integer ngridq(3)
! total number of q-points
integer nqpt
! type of reduction to perform on q-point set (see reducek)
integer reduceq
! number of point group symmetries used for q-point reduction
integer nsymqpt
! point group symmetry matrices used for q-point reduction
integer symqpt(3,3,48)
! locations of q-points on integer grid
integer, allocatable :: ivq(:,:)
! map from non-reduced grid to reduced set
integer, allocatable :: iqmap(:,:,:)
! q-points in lattice coordinates
real(8), allocatable :: vql(:,:)
! q-points in Cartesian coordinates
real(8), allocatable :: vqc(:,:)
! q-point weights
real(8), allocatable :: wqpt(:)
! weights associated with the integral of 1/q^2
real(8), allocatable :: wiq2(:)

end module