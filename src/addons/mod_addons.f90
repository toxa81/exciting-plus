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
complex(8), allocatable :: zdsht(:,:)
! real to complex spherical harmonic transformation
complex(8), allocatable :: yrlm(:,:)
complex(8), allocatable :: dzsht(:,:)
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
! r-vectors of fft grid
real(8), allocatable :: vgrc(:,:)

integer, allocatable :: spnl(:,:)

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
! maximum number of dimensions for the grid
integer ,parameter :: mpigrid_maxndim=3
! actual number of dimensions of the grid
integer mpigrid_ndim
! dimensions of mpi grid
integer mpigrid(mpigrid_maxndim)
! local fraction of fv eigen vectors
complex(8), allocatable :: evecfvloc(:,:,:,:)
! local fraction of sv eigen vectors
complex(8), allocatable :: evecsvloc(:,:,:)
! second-variational Hamiltonian
!complex(8), allocatable :: hmltsvloc(:,:,:)
! dimension for coarse k-point parallelization
integer, parameter :: dim_k=1
integer, parameter :: dim1=1
integer, parameter :: dim2=2
integer, parameter :: dim3=3

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
integer, parameter :: t_seceqnsv_setup=31
integer, parameter :: t_seceqnsv_diag=32

integer, parameter :: t_sic_hunif=35
integer, parameter :: t_sic_genfvprj=36
integer, parameter :: t_sic_wan=37
integer, parameter :: t_sic_wan_gen=38
integer, parameter :: t_sic_wan_ovl=39
integer, parameter :: t_sic_wan_pot=34
integer, parameter :: t_sic_wan_rms=46
integer, parameter :: t_sic_me=47
integer, parameter :: t_sic_wvprod=48

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
