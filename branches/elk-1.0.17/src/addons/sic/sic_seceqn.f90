subroutine sic_seceqn(n_,v_sic,nwork,work_sic,dv_sic)
use modmain
use mod_lf
use mod_nrkp
implicit none
! arguments
integer, intent(in) :: n_
real(8), intent(inout) :: v_sic(n_)
integer, intent(in) :: nwork
real(8), intent(inout) :: work_sic(nwork)
real(8), intent(out) :: dv_sic
! local variables
integer ispn,ist,n,j,ikloc,idm
complex(8), allocatable :: wann_ufv(:,:,:)
logical exist

inquire(file="sic.hdf5",exist=exist)
if (.not.exist) return

allocate(wann_ufv(nwann,nstfv,nspinor))
evalsv=0.d0
do ikloc=1,nkptloc
  wann_ufv=zzero
  do ispn=1,nspinor
    do ist=1,nstfv
      do n=1,nwann
        do j=1,nstsv
          wann_ufv(n,ist,ispn)=wann_ufv(n,ist,ispn)+&
            wann_c(n,j,ikloc)*evecsvloc(ist+(ispn-1)*nstfv,j,ikloc)
        enddo
      enddo
    enddo
  enddo
  call sic_hunif(ikloc,wann_ufv,evecfvloc(1,1,1,ikloc),evecsvloc(1,1,ikloc))
enddo
deallocate(wann_ufv)
call mpi_grid_reduce(evalsv(1,1),nstsv*nkpt,dims=(/dim_k/),all=.true.)
if (wproc) then
! find the occupation numbers and Fermi energy
  call occupy
  if (autoswidth) then
    write(60,*)
    write(60,'("New smearing width : ",G18.10)') swidth
  end if    
! write out the eigenvalues and occupation numbers
  call writeeval
! write the Fermi energy to file
  call writefermi
endif
call mpi_grid_bcast(swidth,dims=(/dim_k,dim2/))
call mpi_grid_bcast(occsv(1,1),nstsv*nkpt,dims=(/dim_k,dim2/))
! set the charge density and magnetisation to zero
rhomt(:,:,:)=0.d0
rhoir(:)=0.d0
if (spinpol) then
  magmt(:,:,:,:)=0.d0
  magir(:,:)=0.d0
end if
do ikloc=1,nkptloc
  call rhomagk(ikloc,evecfvloc(1,1,1,ikloc),evecsvloc(1,1,ikloc))
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
! convert muffin-tin density/magnetisation to spherical harmonics
  call rhomagsh
! symmetrise the density
  call symrf(lradstp,rhomt,rhoir)
! symmetrise the magnetisation
  if (spinpol) call symrvf(lradstp,magmt,magir)
! convert the density from a coarse to a fine radial mesh
  call rfmtctof(rhomt)
! convert the magnetisation from a coarse to a fine radial mesh
  do idm=1,ndmag
    call rfmtctof(magmt(:,:,:,idm))
  end do
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
! compute the effective potential
call poteff
! pack interstitial and muffin-tin effective potential and field into one array
call mixpack(.true.,n_,v_sic)
! mix in the old potential and field with the new
call mixerifc(mixtype,n_,v_sic,dv_sic,nwork,work_sic)
! unpack potential and field
call mixpack(.false.,n_,v_sic)
! add the fixed spin moment effect field
if (fixspin.ne.0) call fsmfield
! Fourier transform effective potential to G-space
call genveffig
return
end