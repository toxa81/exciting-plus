module mod_linresp
implicit none

!-------------------------!
!     Linear response     !
!-------------------------!
! dimension for q-vectors
integer, parameter :: dim_q=2
! dimension for interband transitions
integer, parameter :: dim_b=3
! number of G-vectors for matrix elements
integer ngvecme
! type of linear response calculation
!   0 : charge response
!   1 : magnetic response
integer lrtype
data lrtype/0/

! number of matrix elements <nk|e^{-i(G+q)x}|n'k+q> in the Bloch basis
!  for a given k-point
integer, allocatable :: nmegqblh(:)
! local number of interband transitions (each processor along dim_b does 
!  it's local fraction of transitions)
integer, allocatable :: nmegqblhloc(:,:)
! maximum number of matrix elements <nk|e^{-i(G+q)x}|n'k+q> over all k-points
integer nmegqblhmax
integer nmegqblhlocmax
! matrix elements <nk|e^{-i(G+q)x}|n'k+q> in the Bloch basis
!   1-st index : G-vector
!   2-nd index : global index of pair of bands (n,n')
!   3-rd index : k-point
complex(8), allocatable :: megqblh(:,:,:)
! matrix elements <nk|e^{-i(G+q)x}|n'k+q> in the Bloch basis
!   1-st index : global index of pair of bands (n,n')
!   2-nd index : G-vector
!   3-rd index : k-point
complex(8), allocatable :: megqblh2(:,:)
! pair of bands (n,n') for matrix elements <nk|e^{-i(G+q)x}|n'k+q> by global index
!   1-st index :  1 -> n
!                 2 -> n'
!   2-nd index : global index of pair of bands (n,n')
!   3-rd index : k-point
integer, allocatable :: bmegqblh(:,:,:)

real(8) chi0_include_bands(2)
data chi0_include_bands/-100.1d0,100.1d0/
real(8) chi0_exclude_bands(2)
data chi0_include_bands/100.1d0,-100.1d0/

real(8) lr_min_e12

real(8) megqwan_cutoff1
data megqwan_cutoff1/-0.1d0/
real(8) megqwan_cutoff2
data megqwan_cutoff2/100.1d0/

real(8) megqwan_mindist
data megqwan_mindist/-0.1d0/
real(8) megqwan_maxdist
data megqwan_maxdist/0.1d0/

integer nmegqwan
integer megqwan_tlim(2,3)
integer ntmegqwan
integer, allocatable :: imegqwan(:,:)
integer, allocatable :: itmegqwan(:,:)
integer, allocatable :: idxmegqwan(:,:,:,:,:)
complex(8), allocatable :: megqwan(:,:)
logical :: all_wan_ibt
data all_wan_ibt/.false./

integer nwann_include
data nwann_include/0/
integer, allocatable :: iwann_include(:)

integer nmegqblhwanmax
integer, allocatable :: nmegqblhwan(:)
integer, allocatable :: imegqblhwan(:,:)

complex(8), allocatable :: wann_c_jk(:,:,:)
complex(8), allocatable :: wann_cc(:,:,:)
complex(8), allocatable :: wann_cc2(:,:)

integer ngntujumax
integer, allocatable :: ngntuju(:,:)
integer(2), allocatable :: igntuju(:,:,:,:)
complex(8), allocatable :: gntuju(:,:,:)

complex(8), allocatable :: chi0loc(:,:,:)
complex(8), allocatable :: chi0wanloc(:,:,:)

! array for k+q points
!  1-st index: index of k-point in BZ
!  2-nd index: 1: index of k'=k+q-K
!              2: index of K-vector which brings k+q to first BZ
integer, allocatable :: idxkq(:,:)
! number of energy-mesh points
integer lr_nw
data lr_nw/201/
! first energy point (eV)
real(8) lr_w0
data lr_w0/0.d0/
! last energy point (eV)
real(8) lr_w1
data lr_w1/20.d0/
! energy step
real(8) lr_dw
! energy mesh
complex(8), allocatable :: lr_w(:)
! broadening parameter (eV)
real(8) lr_eta
data lr_eta/0.3d0/

real(8) fxca0
data fxca0/0.d0/
real(8) fxca1
data fxca1/0.d0/
integer nfxca
data nfxca/1/
integer fxctype
data fxctype/0/

! high-level switch: compute chi0 and chi in Wannier functions basis
logical wannier_chi0_chi 
data wannier_chi0_chi/.false./
! high-level switch: .true. if chi0 should be multiplied by 2
logical wannier_chi0_afm
data wannier_chi0_afm/.false./

! low level switch: compute matrix elements of e^{i(G+q)x} in the basis of
!   Wannier functions; depends on crpa and wannier_chi0_chi
logical wannier_megq

! indices of response functions in global array f_response(:,:,:)
integer, parameter :: f_chi0                 = 1
integer, parameter :: f_chi                  = 2
integer, parameter :: f_chi_scalar           = 3
integer, parameter :: f_chi_pseudo_scalar    = 4
integer, parameter :: f_epsilon_matrix_GqGq  = 5
integer, parameter :: f_epsilon_scalar_GqGq  = 6
integer, parameter :: f_inv_epsilon_inv_GqGq = 7
integer, parameter :: f_epsilon_eff          = 8
integer, parameter :: f_epsilon_eff_scalar   = 9
integer, parameter :: f_sigma                = 10
integer, parameter :: f_sigma_scalar         = 11
integer, parameter :: f_loss                 = 12
integer, parameter :: f_loss_scalar          = 13
integer, parameter :: f_chi0_wann            = 14
integer, parameter :: f_chi_wann             = 15
integer, parameter :: f_epsilon_eff_wann     = 16
integer, parameter :: f_sigma_wann           = 17
integer, parameter :: f_loss_wann            = 18

integer, parameter :: nf_response            = 18
complex(8), allocatable :: f_response(:,:,:)

complex(8), allocatable :: uscrnwan(:,:)
complex(8), allocatable :: jscrnwan(:,:)
complex(8), allocatable :: ubarewan(:)
complex(8), allocatable :: u4(:,:,:,:)



end module