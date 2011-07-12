subroutine putevecfd(ikloc,evecfd)
use modmain
use mod_mpi_grid
implicit none
integer, intent(in) :: ikloc
complex(8), intent(in) :: evecfd(nspinor*nmatmax,nstsv)
!
integer ik,n
!
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
inquire(iolength=n) vkl(:,ik),vgkl(:,:,1,ikloc),igkig(:,1,ikloc),&
  ngkmax,nmatmax,nspinor,nstsv,evecfd
open(70,file=trim(scrpath)//'EVECFD'//trim(filext),action='WRITE', &
  form='UNFORMATTED',access='DIRECT',recl=n)
write(70,rec=ik) vkl(:,ik),vgkl(:,:,1,ikloc),igkig(:,1,ikloc),&
  ngkmax,nmatmax,nspinor,nstsv,evecfd
close(70)
return
end subroutine
