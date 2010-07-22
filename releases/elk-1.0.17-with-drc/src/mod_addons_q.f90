module mod_addons_q
!-------------------------!
!     q and G+q vectors   !
!-------------------------!
! number of q-vectors
integer nvq
! q-vectors in k-mesh coordinates
integer, allocatable :: vqm(:,:)
! non-reduced (to first BZ) q-vectors in lattice coordinates
real(8), allocatable :: vqlnr(:,:)
! non-reduced (to first BZ) q-vectors in Cartesian coordinates
real(8), allocatable :: vqcnr(:,:)
! reduced to first BZ q-vectors in lattice coordinates
real(8), allocatable :: vql(:,:)
! reduced to first BZ q-vectors in Cartesian coordinates
real(8), allocatable :: vqc(:,:)
! cutoff for |G+q|
real(8) gqmax
data gqmax/2.d0/
integer gqsh
data gqsh/2/
! maximum number of G+q vectors
integer ngqmax
! index of G0-vector, which brigs q-vector to first BZ
integer, allocatable :: ig0q(:)
integer iig0q
! number of G+q vectors
integer, allocatable :: ngq(:)
! G+q vectors in Cartesian coordinates
real(8), allocatable ::  vgqc(:,:,:)
! length of |G+q| vectors
real(8), allocatable :: gq(:,:)
! global index of G of G+q vector  
integer, allocatable :: igqig(:,:)
! 4*Pi/|G+q|^2 (Fourier transform of Hartree potential)
real(8), allocatable :: vhgq(:,:)

real(8), allocatable :: tpgq(:,:)
complex(8), allocatable :: sfacgq(:,:)
complex(8), allocatable :: ylmgq(:,:)

! true if Gamma point included in the integration over Brillouin zone
logical tq0bz
data tq0bz/.true./
integer nvq0
real(8) vq0c(3,8)
real(8) aq0(8)


end module