subroutine poteff1
use modmain
implicit none
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
integer i,ik,ikloc
real(8) spzn1(maxspecies)
integer np
real(8), allocatable :: vpl(:,:)
real(8), allocatable :: fp(:)

occsv=0.d0
! three p-bands of MgO 
occsv(6:8,:)=1.d0
rhomt(:,:,:)=0.d0
rhoir(:)=0.d0
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
! read eigen-vectors
if (mpi_grid_side(dims=(/dim_k/))) then
  do i=0,mpi_grid_size(dim_k)-1
    if (i.eq.mpi_grid_x(dim_k)) then
      do ikloc=1,nkptloc
        ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
        call getevecfv(vkl(1,ik),vgkl(1,1,1,ikloc),evecfv)
        call getevecsv(vkl(1,ik),evecsv)
        call rhomagk(ikloc,evecfv,evecsv)
      enddo !ikloc
    endif
    if (.not.parallel_read) call mpi_grid_barrier(dims=(/dim_k/))
  enddo
endif !mpi_grid_side(dims=(/dim_k/)
call mpi_grid_reduce(rhomt(1,1,1),lmmaxvr*nrmtmax*natmtot,dims=(/dim_k,2/),&
  all=.true.)
call mpi_grid_reduce(rhoir(1),ngrtot,dims=(/dim_k,2/),all=.true.)
call rhomagsh
call symrf(lradstp,rhomt,rhoir)
call rfmtctof(rhomt)
spzn1=spzn
spzn=0.d0
call poteff
spzn=spzn1 

np=500
allocate(vpl(3,np))
do i=1,np
  vpl(:,i)=1.d0*(/i,i,i/)/np
enddo
allocate(fp(np))
call rfarray(lmaxvr,lmmaxvr,veffmt,veffir,np,vpl,fp)
do i=1,np
  write(100,*)i,fp(i)
enddo


return
end