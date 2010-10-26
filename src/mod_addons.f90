module mod_addons

! number of atom classes (non-equivalent atoms)
integer natmcls
! i-th class -> ias mapping
integer, allocatable :: ic2ias(:)
! ias -> ic mapping
integer, allocatable :: ias2ic(:)
! ias -> is mapping
integer, allocatable :: ias2is(:)
! ias -> ia mapping
integer, allocatable :: ias2ia(:)
! lm -> l mapping
integer, allocatable :: lm2l(:)
! lm -> m mapping
integer, allocatable :: lm2m(:)
! maximum number of muffin-tin radial functions
integer nufrmax
! muffin-tin radial functions
!  1-st index: radial point (from 1 to nrmtmax)
!  2-nd index: l (from 0 to lmaxvr)
!  3-rd index: order of function (from 1 to nfrmax; apwfr first, lofr then)
!  4-th index: atom class (from 1 to natmcls)
real(8), allocatable :: ufr(:,:,:,:)
! product <u_l^{io1}|u_l^{io2}> of two radial functions in the muffin-tin
!  1-st index: l (from 0 to lmaxvr)
!  2-nd index: order of first function
!  3-rd index: order of second function
!  4-th index: atom class (from 1 to natmcls)
real(8), allocatable :: ufrp(:,:,:,:)
! number of radial functions for a given l
!  1-st index: l (from 0 to lmaxvr)
!  2-nd index: atom species (from 1 to nspecies)
integer, allocatable :: nufr(:,:)
! number of atoms with local point symmetry
integer natlps
! transformation matrices for real spherical harmonics
real(8), allocatable :: lpsrsh(:,:,:)
! index of atom which has local point symmetry
integer, allocatable :: iatlps(:)
! complex to real spherical harmonic transformation
complex(8), allocatable :: rylm(:,:)
! real to complex spherical harmonic transformation
complex(8), allocatable :: yrlm(:,:)
! global complex to real lps spherical harmonic transformation
complex(8), allocatable :: rylm_lps(:,:,:)
! real lps to global complex spherical harmonic transformation
complex(8), allocatable :: yrlm_lps(:,:,:)
! .true. if density matrix is computed instead of occupancy matrix
logical ldensmtrx
! energy (or band) interval for bands for density matrix 
real(8) dm_e1,dm_e2
! number of nearest neighbours for each atom
integer, allocatable :: nnghbr(:)
! list of nearest neighbours
integer, allocatable :: inghbr(:,:,:)
! unit conversion
real(8), parameter :: ha2ev=27.21138386d0
real(8), parameter :: au2ang=0.5291772108d0
! super-cell vectors in lattice coordinates
integer scvl(3,3)

integer debug_level
data debug_level/0/
integer fdbgout
character*256 fdbgname

!-----------------------!
!      MPI parallel     !
!-----------------------!
! local number of k-points
integer nkptloc
! local number of non-reduced k-points
integer nkptnrloc
! .true. if processor writes some info
logical wproc
! .true. if mpi grid layout comes from input file
logical lmpigrid
data lmpigrid/.false./
! dimensions of mpi grid
integer mpigrid(3)
! local fraction of fv eigen vectors
complex(8), allocatable :: evecfvloc(:,:,:,:)
! local fraction of sv eigen vectors
complex(8), allocatable :: evecsvloc(:,:,:)
! dimension for coarse k-point parallelization
integer, parameter :: dim_k=1
! auxiliary dimension for fine parallelization 
integer, parameter :: dim2=2

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

integer nmegqwanmax
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

!------------------!
!     Wannier      !
!------------------!
logical wannier
integer wann_natom
integer wann_norbgrp
integer wann_ntype
logical wann_add_poco
integer, allocatable :: wann_norb(:)
integer, allocatable :: wann_iorb(:,:,:)
integer, allocatable :: wann_iprj(:,:)
real(8), allocatable :: wann_eint(:,:)
real(8), allocatable :: wann_v(:)

integer nwann
integer, allocatable :: iwann(:,:)
integer, allocatable :: nwannias(:)
integer nwannloc
  
! expansion coefficients of Wannier functions over spinor Bloch eigen-functions  
complex(8), allocatable :: wann_c(:,:,:)
! Bloch-sums of WF
complex(8), allocatable :: wann_unkmt(:,:,:,:,:,:)
complex(8), allocatable :: wann_unkit(:,:,:,:)

! H(k) in WF basis
complex(8), allocatable :: wann_h(:,:,:)
! e(k) of WF H(k) (required for band-sctructure plot only)
real(8), allocatable :: wann_e(:,:)
! momentum operator in WF basis
complex(8), allocatable :: wann_p(:,:,:,:)

real(8), allocatable :: wann_ene(:)
real(8), allocatable :: wann_occ(:)

real(8) zero3d(3)
real(8) bound3d(3,3)
integer nrxyz(3)
integer nwfplot
integer firstwf
logical wannier_lc
integer nwann_lc
integer, allocatable :: wann_iorb_lc(:,:,:)
real(8), allocatable :: wann_iorb_lcc(:,:)

integer nwann_h
integer, allocatable :: iwann_h(:)

logical wannier_soft_eint
data wannier_soft_eint/.false./
real(8), allocatable :: wannier_soft_eint_w1(:)
real(8), allocatable :: wannier_soft_eint_w2(:)
real(8), allocatable :: wannier_soft_eint_e1(:)
real(8), allocatable :: wannier_soft_eint_e2(:)
real(8) wannier_min_prjao
data wannier_min_prjao/-0.1d0/

logical ldisentangle
data ldisentangle/.false./

!----------------!
!      timer     !
!----------------!
integer, parameter :: t_iter_tot=2
integer, parameter :: t_init=10

integer, parameter :: t_seceqn=18
integer, parameter :: t_seceqnfv=19
integer, parameter :: t_seceqnfv_setup=20
integer, parameter :: t_seceqnfv_setup_h=21
integer, parameter :: t_seceqnfv_setup_h_mt=22
integer, parameter :: t_seceqnfv_setup_h_it=23
integer, parameter :: t_seceqnfv_setup_o=24
integer, parameter :: t_seceqnfv_setup_o_mt=25
integer, parameter :: t_seceqnfv_setup_o_it=26
integer, parameter :: t_seceqnfv_diag=27

integer, parameter :: t_seceqnsv=30
integer, parameter :: t_svhmlt_setup=31
integer, parameter :: t_svhmlt_diag=32
integer, parameter :: t_svhmlt_tot=33

integer, parameter :: t_apw_rad=40
integer, parameter :: t_rho_mag_sum=41
integer, parameter :: t_rho_mag_sym=42
integer, parameter :: t_rho_mag_tot=43
integer, parameter :: t_pot=44
integer, parameter :: t_dmat=45

!-------------!
!      SIC    !
!-------------!
logical sic
data sic/.false./
! product of a Wannier function with it's potential
complex(8), allocatable :: wvmt(:,:,:,:,:,:)
complex(8), allocatable :: wvir(:,:,:,:)
! Wannier functions
complex(8), allocatable :: wanmt(:,:,:,:,:,:)
complex(8), allocatable :: wanir(:,:,:,:)
! matrix elements of Wannier potential <W_{n0}|V_{n0}|W_{n'T}>
complex(8), allocatable :: vwanme(:)
! LDA Hamiltonian in k-space in Wannier basis 
complex(8), allocatable :: sic_wann_h0k(:,:,:)
! LDA energies of Wannier functions
real(8), allocatable :: sic_wann_e0(:)
integer nsclsic
data nsclsic/3/
integer isclsic
real(8) sic_etot_correction
data sic_etot_correction/0.d0/
real(8) wann_r_cutoff
data wann_r_cutoff/6.0/
complex(8), allocatable :: sic_wb(:,:,:,:)
complex(8), allocatable :: sic_wvb(:,:,:,:)


!--------------!
!      PAPI    !
!--------------!
integer, parameter :: maxpapievents=16
character*256 :: papievent(maxpapievents)
data papievent(1)/"PAPI_FP_OPS"/
integer :: npapievents
data npapievents/1/

integer, parameter :: pt_resp_tot=1
integer, parameter :: pt_megq=2
integer, parameter :: pt_megqblh=3
integer, parameter :: pt_megqblh2=4
integer, parameter :: pt_chi0_zgemm=5
integer, parameter :: pt_chi0=6
integer, parameter :: pt_chi=7
integer, parameter :: pt_crpa_tot1=8
integer, parameter :: pt_crpa_tot2=9
integer, parameter :: pt_uscrn=10
integer, parameter :: pt_vscrn=11
integer, parameter :: pt_megqblh_mt=12
integer, parameter :: pt_megqblh_it=13

end module
