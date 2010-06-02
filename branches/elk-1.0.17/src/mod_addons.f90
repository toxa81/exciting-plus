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
!integer nvq0
!integer, allocatable :: ivq0m_list(:,:)
! q-vector in lattice coordinates
!real(8) vq0l(3)
! q-vector in Cartesian coordinates
!real(8) vq0c(3)
! reduced q-vector in lattice coordinates
!real(8) vq0rl(3)
! reduced q-vector in Cartesian coordinates
!real(8) vq0rc(3)
! index of G-vector which brings q to first BZ
integer lr_igq0
! first G-shell for matrix elements
integer gshme1
! last G-shell for matrix elements
integer gshme2
! first G-vector for matrix elements
integer gvecme1
! last G-vector for matrix elements
integer gvecme2
! number of G-vectors for matrix elements
integer ngvecme
real(8) maxomega
real(8) domega
real(8) lr_eta
real(8) lr_e1,lr_e2
! type of linear response calculation
!   0 : charge response
!   1 : magnetic response
integer lrtype
real(8) lr_min_e12
real(8) lr_e1_wan
real(8) lr_e2_wan

! G+q vectors in Cart.coord.
real(8), allocatable :: lr_vgq0c(:,:)
! length of G+q vectors
real(8), allocatable :: lr_gq0(:)
! theta and phi angles of G+q vectors
real(8), allocatable :: lr_tpgq0(:,:)
! sperical harmonics of G+q vectors
complex(8), allocatable :: lr_ylmgq0(:,:)
! structure factor for G+q vectors
complex(8), allocatable :: lr_sfacgq0(:,:)

! number of matrix elements <nk|e^{-i(G+q)x}|n'k+q> in the Bloch basis
!  for a given k-point
integer, allocatable :: nmegqblh(:)
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

logical megqwan_afm
data megqwan_afm/.false./

integer nmegqwanmax
integer nmegqwan
integer megqwan_tlim(2,3)
integer, allocatable :: imegqwan(:,:)
integer, allocatable :: idxmegqwan(:,:,:,:,:)
complex(8), allocatable :: megqwan(:,:)

integer nmegqblhwanmax
integer, allocatable :: nmegqblhwan(:)
integer, allocatable :: imegqblhwan(:,:)

complex(8), allocatable :: wann_cc(:,:,:)
complex(8), allocatable :: wann_cc2(:,:)


integer ngntujumax
integer, allocatable :: ngntuju(:,:)
integer(2), allocatable :: igntuju(:,:,:,:)
complex(8), allocatable :: gntuju(:,:,:)




! array for k and k+q stuff
!  1-st index: index of k-point in BZ
!  2-nd index: 1: index of k'=k+q-K
!              2: index of K-vector which brings k+q to first BZ
integer, allocatable :: idxkq(:,:)
! number of energy-mesh points
integer nepts
! energy mesh
complex(8), allocatable :: lr_w(:)
real(8) fxca0
real(8) fxca1
integer nfxca
integer fxctype

! high-level switch: solve scalar equation for chi
logical scalar_chi
data scalar_chi/.false./
! high-level switch: split file with matrix elements over k-points
logical split_megq_file
data split_megq_file/.false./
! high-level switch:: read files in parallel
logical parallel_read
data parallel_read/.true./
! high-level switch:: write files in parallel (where it is possible)
logical parallel_write
data parallel_write/.true./
! high-level switch: compute chi0 and chi in Wannier functions basis
logical wannier_chi0_chi 
data wannier_chi0_chi/.false./
! low level switch: compute matrix elements of e^{i(G+q)x} in the basis of
!   Wannier functions; depends on crpa and wannier_chi0_chi
logical wannier_megq
! low-level switch: write or not file with matrix elements; depends on task 
logical write_megq_file
! low level switch: compute screened W matrix; depends on crpa
logical screened_w
data screened_w/.false./
logical screened_u
data screened_u/.false./
logical write_chi0_file

real(8) megqwan_maxdist

logical crpa
real(8) crpa_e1,crpa_e2

integer, allocatable :: spinor_ud(:,:,:)

real(8), allocatable :: lr_occsvnr(:,:)
real(8), allocatable :: lr_evalsvnr(:,:)

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
integer, parameter :: f_chi0_wann_full       = 14
integer, parameter :: f_chi0_wann            = 15
integer, parameter :: f_chi_wann             = 16
integer, parameter :: f_epsilon_eff_wann     = 17
integer, parameter :: f_sigma_wann           = 18
integer, parameter :: f_loss_wann            = 19

integer, parameter :: nf_response            = 19
complex(8), allocatable :: f_response(:,:,:)

complex(8), allocatable :: uscrnwan(:,:,:)
complex(8), allocatable :: ubarewan(:,:)









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

complex(8), allocatable :: wf_v_mtrx(:,:,:,:,:)

real(8) zero3d(3)
real(8) bound3d(3,3)
integer nrxyz(3)
integer nwfplot
integer firstwf
integer iwfv
logical wannier_lc
integer nwann_lc
integer, allocatable :: wann_iorb_lc(:,:,:)
real(8), allocatable :: wann_iorb_lcc(:,:)

integer nwann_h
integer, allocatable :: iwann_h(:)

logical wannier_soft_eint
real(8) wannier_soft_eint_width
real(8) wannier_soft_eint_e1
real(8) wannier_soft_eint_e2
real(8) wannier_min_prjao

logical ldisentangle

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

end module