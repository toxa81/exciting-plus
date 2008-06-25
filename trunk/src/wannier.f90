module   wannier
implicit none

  integer                :: wf_dim
  integer   ,allocatable :: wf_n(:,:)

contains

  subroutine wf_init
  implicit   none
  
  wf_dim = 3
  allocate(wf_n(wf_dim,6))
    
  
  end subroutine wf_init 













end module
