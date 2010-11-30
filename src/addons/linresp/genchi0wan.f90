subroutine genchi0wan(mexp,chi0wan_k,chi0wan)
use modmain
use mod_linresp
use mod_wannier
implicit none
complex(8), intent(in) :: chi0wan_k(nmegqwan,nmegqwan,nkptnrloc)
complex(8), intent(in) :: mexp(nmegqwan,nmegqwan,nkptnrloc) !nkptnrloc,ntrchi0wan)
complex(8), intent(out) :: chi0wan(nmegqwan,nmegqwan)
integer ikloc

! TODO: utilize second dimension
chi0wan(:,:)=zzero
do ikloc=1,nkptnrloc
  chi0wan(:,:)=chi0wan(:,:)+mexp(:,:,ikloc)*chi0wan_k(:,:,ikloc)
enddo
! sum chi0wan over all k-points
call mpi_grid_reduce(chi0wan(1,1),nmegqwan*nmegqwan,dims=(/dim_k/))
chi0wan(:,:)=chi0wan(:,:)/nkptnr/omega
return
end