subroutine genwann(ikloc,evecfv,evecsv)
use modmain
use mod_mpi_grid
implicit none
! arguments
integer, intent(in) :: ikloc
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
! local variables
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfsvmt(:,:,:,:,:)
complex(8), allocatable :: wfsvit(:,:,:)
integer :: ik

ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
! allocate arrays
allocate(wfsvmt(lmmaxvr,nufrmax,natmtot,nspinor,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
call match(ngk(1,ik),gkc(1,1,ikloc),tpgkc(1,1,1,ikloc),sfacgk(1,1,1,ikloc),apwalm)
! generate second-varioational wave-functions
call genwfsvmt(lmaxvr,lmmaxvr,ngk(1,ik),evecfv,evecsv,apwalm,wfsvmt)
! calculate WF expansion coefficients
call genwann_c(ik,vkc(:,ik),evalsv(1,ik),wfsvmt,wann_c(1,1,ikloc))
deallocate(wfsvmt,apwalm)
return
end


