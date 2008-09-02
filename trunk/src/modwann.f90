module modwann
implicit none

logical wannier

integer wf_dim
integer, allocatable :: wf_n(:,:)
  
integer, allocatable :: wf_lhbnd(:,:)
  
complex(8), allocatable :: a_ort(:,:,:,:)
complex(8), allocatable :: wf_h(:,:,:,:)
real(8), allocatable :: wf_e(:,:,:)
  
end module
