module modwann
implicit none

integer wann_natoms
integer wann_nspins
logical wann_use_lhen
logical wann_add_poco
integer, allocatable :: wann_iatom(:)
integer, allocatable :: wann_iorb(:,:)
real(8), allocatable :: wann_lhen(:,:,:,:)
integer, allocatable :: wann_lhbnd(:,:,:,:)
real(8), allocatable :: wann_deltav(:,:,:)

integer wf_dim
integer, allocatable :: wf_n(:,:)
  
integer, allocatable :: wf_lhbnd(:,:,:)
real(8), allocatable :: wf_lhen(:,:,:)
real(8) wf_e1,wf_e2
integer wf_n1,wf_n2
  
!complex(8), allocatable :: a_ort(:,:,:,:)
complex(8), allocatable :: wfc(:,:,:,:)
complex(8), allocatable :: wf_h(:,:,:,:)
real(8), allocatable :: wf_e(:,:,:)
real(8), allocatable :: wf_deltav(:,:)

complex(8), allocatable :: wfpoco(:,:,:)
complex(8), allocatable :: wfpoco1(:,:,:)

complex(8), allocatable :: psao(:)
  
end module
