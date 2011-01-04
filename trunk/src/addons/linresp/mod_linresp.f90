module mod_linresp
implicit none

! dimension for q-vectors
integer, parameter :: dim_q=2

! type of linear response calculation
!   0 : charge response
!   1 : magnetic response (not yet implemented)
integer lrtype
data lrtype/0/

! temporary array used in genchi0blh
complex(8), allocatable :: megqblh2(:,:)

! interval of bands to include in chi0
real(8) chi0_include_bands(2)
data chi0_include_bands/-100.1d0,100.1d0/

! interval of bands to exclude from chi0
real(8) chi0_exclude_bands(2)
data chi0_exclude_bands/100.1d0,-100.1d0/

complex(8), allocatable :: wann_cc(:,:,:)
complex(8), allocatable :: wann_cc2(:,:)

complex(8), allocatable :: chi0loc(:,:,:)
complex(8), allocatable :: chi0wanloc(:,:,:)

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

complex(8), allocatable :: u4(:,:,:,:)
logical screenu4
data screenu4/.true./

end module
