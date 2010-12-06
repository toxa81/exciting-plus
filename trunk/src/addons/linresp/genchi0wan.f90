subroutine genchi0wan(mexp,chi0wan_k,chi0wan)
use modmain
use mod_expigqr
implicit none
complex(8), intent(in) :: chi0wan_k(megqwantran%nwt,megqwantran%nwt,nkptnrloc)
complex(8), intent(in) :: mexp(megqwantran%nwt,megqwantran%nwt,nkptnrloc)
complex(8), intent(out) :: chi0wan(megqwantran%nwt,megqwantran%nwt)
integer ikloc

! TODO: utilize second dimension
chi0wan(:,:)=zzero
do ikloc=1,nkptnrloc
  chi0wan(:,:)=chi0wan(:,:)+mexp(:,:,ikloc)*chi0wan_k(:,:,ikloc)
enddo
! sum chi0wan over all k-points
call mpi_grid_reduce(chi0wan(1,1),megqwantran%nwt*megqwantran%nwt,dims=(/dim_k/))
chi0wan(:,:)=chi0wan(:,:)/nkptnr/omega
return
end