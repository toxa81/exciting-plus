! full diagonalization
subroutine seceqnfd(ikloc,evecfd)
use modmain
use mod_seceqn
use mod_sic
implicit none
integer, intent(in) :: ikloc
complex(8), intent(out) :: evecfd(nspinor*nmatmax,nstsv)
integer ik,i,ispn,j
complex(8) zt1
! allocatable arrays
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: h(:,:,:)
complex(8), allocatable :: o(:,:,:)
!
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
call timer_start(t_seceqn)
call timer_start(t_seceqnfv)
call timer_start(t_seceqnfv_setup)
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
! find the matching coefficients
call match(ngk(1,ik),gkc(:,1,ikloc),tpgkc(:,:,1,ikloc),&
  sfacgk(:,:,1,ikloc),apwalm)
allocate(h(nmat(1,ik),nmat(1,ik),nspinor))
allocate(o(nmat(1,ik),nmat(1,ik),nspinor))
h=zzero
o=zzero
call setovl(ngk(1,ik),nmat(1,ik),igkig(1,1,ikloc),apwalm,o)
if (spinpol) then
  o(:,:,2)=o(:,:,1)
endif
if (spinpol) then
  do ispn=1,nspinor
    call sethml(ngk(1,ik),nmat(1,ik),vgkc(1,1,1,ikloc),igkig(1,1,ikloc),&
      apwalm,h(1,1,ispn),ispn,ispn)
  enddo
else
  call sethml(ngk(1,ik),nmat(1,ik),vgkc(1,1,1,ikloc),igkig(1,1,ikloc),&
    apwalm,h)
endif
if (sic.and..not.tsicsv) then
  call sic_genbprj(ikloc,apwalm=apwalm)
  call sic_hunif_fd(ikloc,nmat(1,ik),h,o)
endif
deallocate(apwalm)
call timer_stop(t_seceqnfv_setup)
call timer_start(t_seceqnfv_diag)
evecfd=zzero
if (mpi_grid_root((/dim2/))) then
  do ispn=1,nspinor
    i=(ispn-1)*nstfv+1
    j=(ispn-1)*nmatmax+1
    call diagzheg(nmat(1,ik),nstfv,nspinor*nmatmax,evaltol,&
      h(1,1,ispn),o(1,1,ispn),evalsv(i,ik),evecfd(j,i))
  enddo
endif
call timer_stop(t_seceqnfv_diag)
deallocate(o,h)
call timer_start(t_seceqnfv)
call timer_stop(t_seceqn)
if (sic.and.tsicsv) then
  call sic_hunif_sv(ikloc,evecfd)
endif
if (wannier) then
  call wan_gencsv_aux(ikloc,evecfd=evecfd)
endif
return
end subroutine

