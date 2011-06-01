subroutine rhomag1
use modmain
implicit none
integer ikloc,idm
call timer_start(t_rho_mag_tot)
! set the charge density and magnetisation to zero
rhomt(:,:,:)=0.d0
rhoir(:)=0.d0
if (spinpol) then
  magmt(:,:,:,:)=0.d0
  magir(:,:)=0.d0
end if
do ikloc=1,nkptloc
  call rhomagk1(ikloc,evecfvloc(1,1,1,ikloc),evecsvloc(1,1,ikloc))
end do
call mpi_grid_reduce(rhomt(1,1,1),lmmaxvr*nrmtmax*natmtot,&
  dims=(/dim_k,dim2/))
call mpi_grid_reduce(rhoir(1),ngrtot,dims=(/dim_k,dim2/))
if (spinpol) then
  call mpi_grid_reduce(magmt(1,1,1,1),lmmaxvr*nrmtmax*natmtot*ndmag,&
    dims=(/dim_k,dim2/))
  call mpi_grid_reduce(magir(1,1),ngrtot*ndmag,dims=(/dim_k,dim2/))
endif
if (mpi_grid_root(dims=(/dim_k,dim2/))) then
! symmetrise the density
  call symrf(1,rhomt,rhoir)
! symmetrise the magnetisation
  if (spinpol) call symrvf(1,magmt,magir)
! add the core density to the total density
  call addrhocr
endif
call mpi_grid_bcast(rhomt(1,1,1),lmmaxvr*nrmtmax*natmtot,&
  dims=(/dim_k,dim2/))
call mpi_grid_bcast(rhoir(1),ngrtot,dims=(/dim_k,dim2/))
if (spinpol) then
  call mpi_grid_bcast(magmt(1,1,1,1),lmmaxvr*nrmtmax*natmtot*ndmag,&
    dims=(/dim_k,dim2/))
  call mpi_grid_bcast(magir(1,1),ngrtot*ndmag,dims=(/dim_k,dim2/))
endif
! calculate the charges
call charge
! calculate the moments
if (spinpol) call moment
! normalise the density
call rhonorm
call timer_stop(t_rho_mag_tot)
return
end subroutine
