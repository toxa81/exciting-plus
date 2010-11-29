module mod_wannier
implicit none

type wannier_transition
! number of Wannier functions taken
  integer nwan
! mapping to global index
  integer, allocatable :: idxwan(:)
  integer nwantran
  integer, allocatable :: iwantran(:,:) 
end type wannier_transition



end module
