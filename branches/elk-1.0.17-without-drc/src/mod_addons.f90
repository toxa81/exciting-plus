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

end module