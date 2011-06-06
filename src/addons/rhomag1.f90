subroutine rhomag1
use modmain
implicit none
integer ikloc,idm
integer ias,ic,is,lm3,l1,io1,j1,l2,io2,j2,ispn
real(8), allocatable :: f1(:,:)
call timer_start(t_rho_mag_tot)
! set the charge density and magnetisation to zero
rhomt(:,:,:)=0.d0
rhoir(:)=0.d0
if (spinpol) then
  magmt(:,:,:,:)=0.d0
  magir(:,:)=0.d0
end if
rhomagmt=0.d0
do ikloc=1,nkptloc
  call rhomagk1(ikloc,evecfvloc(1,1,1,ikloc),evecsvloc(1,1,ikloc))
end do
call mpi_grid_reduce(rhomagmt(1,1,1,1,1),&
  nlufrmax*nlufrmax*lmmaxvr*natmtot*nspinor,dims=(/dim_k/))
allocate(f1(nrmtmax,nspinor))
do ias=1,natmtot
  ic=ias2ic(ias)
  is=ias2is(ias)
  do lm3=1,lmmaxvr
    f1=0.d0
    j1=0
    do l1=0,lmaxvr
      do io1=1,nufr(l1,is)
        j1=j1+1
        j2=0
        do l2=0,lmaxvr
          do io2=1,nufr(l2,is)
            j2=j2+1
            do ispn=1,nspinor
              f1(:,ispn)=f1(:,ispn)+rhomagmt(j1,j2,lm3,ias,ispn)*ufr(:,l1,io1,ic)*&
                ufr(:,l2,io2,ic)
            enddo
          enddo
        enddo !l2
      enddo 
    enddo !l1
    if (spinpol) then
      rhomt(lm3,:,ias)=f1(:,1)+f1(:,2)
      magmt(lm3,:,ias,1)=f1(:,1)-f1(:,2)
    endif
  enddo !lm3
enddo !ias
deallocate(f1)
!call mpi_grid_reduce(rhomt(1,1,1),lmmaxvr*nrmtmax*natmtot,&
!  dims=(/dim_k,dim2/))
call mpi_grid_reduce(rhoir(1),ngrtot,dims=(/dim_k,dim2/))
if (spinpol) then
  !call mpi_grid_reduce(magmt(1,1,1,1),lmmaxvr*nrmtmax*natmtot*ndmag,&
  !  dims=(/dim_k,dim2/))
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
