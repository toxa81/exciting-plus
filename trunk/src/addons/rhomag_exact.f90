subroutine rhomag_exact
use modmain
implicit none
integer ikloc,idm,n
integer ias,ic,is,lm3,l1,io1,j1,l2,io2,j2,ispn
real(8), allocatable :: fr(:,:)
call timer_start(t_rho_mag_tot)
! set the charge density and magnetisation to zero
rhomt(:,:,:)=0.d0
rhoir(:)=0.d0
if (spinpol) then
  magmt(:,:,:,:)=0.d0
  magir(:,:)=0.d0
end if
rhomagmt=0.d0
rhomagit=zzero
call timer_start(t_rho_mag_sum)
do ikloc=1,nkptloc
  call rhomagk_exact(ikloc,evecfvloc(1,1,1,ikloc),evecsvloc(1,1,ikloc))
end do
n=nlufrmax*nlufrmax*lmmaxvr*natmtot*nspinor
call mpi_grid_reduce(rhomagmt(1,1,1,1,1),n,dims=(/dim_k/))
call mpi_grid_reduce(rhomagit(1,1),ngrtot*nspinor,dims=(/dim_k/))
call timer_stop(t_rho_mag_sum)
call timer_start(t_rho_mag_conv)
allocate(fr(nrmtmax,nspinor))
do ias=1,natmtot
  ic=ias2ic(ias)
  is=ias2is(ias)
  do lm3=1,lmmaxvr
    fr=0.d0
    j1=0
    do l1=0,lmaxvr
      do io1=1,nufr(l1,is)
        j1=j1+1
        j2=0
        do l2=0,lmaxvr
          do io2=1,nufr(l2,is)
            j2=j2+1
            do ispn=1,nspinor
              fr(:,ispn)=fr(:,ispn)+rhomagmt(j1,j2,lm3,ias,ispn)*ufr(:,l1,io1,ic)*&
                ufr(:,l2,io2,ic)
            enddo
          enddo
        enddo !l2
      enddo 
    enddo !l1
    if (spinpol) then
      rhomt(lm3,:,ias)=fr(:,1)+fr(:,2)
      magmt(lm3,:,ias,1)=fr(:,1)-fr(:,2)
    else
      rhomt(lm3,:,ias)=fr(:,1)
    endif
  enddo !lm3
enddo !ias
deallocate(fr)
do ispn=1,nspinor
  call zfftifc(3,ngrid,1,rhomagit(1,ispn))
enddo
if (spinpol) then
  rhoir(:)=rhoir(:)+dreal(rhomagit(:,1)+rhomagit(:,2))
  magir(:,1)=magir(:,1)+dreal(rhomagit(:,1)-rhomagit(:,2))
else
  rhoir(:)=rhoir(:)+dreal(rhomagit(:,1))
endif
call timer_stop(t_rho_mag_conv)
if (mpi_grid_root(dims=(/dim_k,dim2/))) then
  call timer_start(t_rho_mag_sym)
! symmetrise the density
  call symrf(1,rhomt,rhoir)
! symmetrise the magnetisation
  if (spinpol) call symrvf(1,magmt,magir)
  call timer_stop(t_rho_mag_sym)
! add the core density to the total density
  call addrhocr
endif
call mpi_grid_bcast(rhomt(1,1,1),lmmaxvr*nrmtmax*natmtot)
call mpi_grid_bcast(rhoir(1),ngrtot)
if (spinpol) then
  call mpi_grid_bcast(magmt(1,1,1,1),lmmaxvr*nrmtmax*natmtot*ndmag)
  call mpi_grid_bcast(magir(1,1),ngrtot*ndmag)
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
