subroutine sic
use modmain
implicit none

integer i1,i2,i3,n
integer ik,ikloc,ik1,j,ig,sz,i,isym
integer n1,n2
integer ntr,itr,ntrloc,itrloc
integer, allocatable :: vtrl(:,:)
integer, allocatable :: ivtit(:,:,:)
integer, parameter :: trmax=1

complex(8), allocatable :: wanmt(:,:,:,:,:,:)
complex(8), allocatable :: wanir(:,:,:,:)
complex(8), allocatable :: wanmt0(:,:,:,:,:)
complex(8), allocatable :: wanir0(:,:,:)

integer, allocatable :: ngknr(:)
integer, allocatable :: igkignr(:,:)
real(8), allocatable :: vgklnr(:,:,:)
real(8), allocatable :: vgkcnr(:,:,:)
real(8), allocatable :: gknr(:,:)
real(8), allocatable :: tpgknr(:,:,:)
complex(8), allocatable :: sfacgknr(:,:,:)
complex(8), allocatable :: ylmgknr(:,:,:)

complex(8), allocatable :: wfsvmtloc(:,:,:,:,:,:)
complex(8), allocatable :: wfsvitloc(:,:,:,:)
complex(8), allocatable :: wfsvcgloc(:,:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: ovlm(:,:)
complex(8), external :: zfinp_
complex(8), external :: inner_product


! mpi grid layout
!          (2)
!     +----+----+--> T-vectos 
!     |    |    |
!     +----+----+--
! (1) |    |    |
!     +----+----+--
!     |    |    |
!     v
!  k-points


! initialise universal variables
call init0
call init1
if (.not.mpi_grid_in()) return
wproc=mpi_grid_root()
! read the density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
call geturf
call genurfprod
! read Fermi energy
if (mpi_grid_root()) call readfermi
call mpi_grid_bcast(efermi)

if (wproc) then
  open(151,file='SIC.OUT',form='FORMATTED',status='REPLACE')
endif

! get energies of states in reduced part of BZ
call timer_start(3,reset=.true.)
if (wproc) then
  write(151,*)
  write(151,'("Reading energies of states")')
  call flushifc(151)
! read from IBZ
  do ik=1,nkpt
    call getevalsv(vkl(1,ik),evalsv(1,ik))
  enddo
endif
call mpi_grid_bcast(evalsv(1,1),nstsv*nkpt)
allocate(lr_evalsvnr(nstsv,nkptnr))
lr_evalsvnr=0.d0
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  call findkpt(vklnr(1,ik),isym,ik1) 
  lr_evalsvnr(:,ik)=evalsv(:,ik1)
enddo
call timer_stop(3)
if (wproc) then
  write(151,'("Done in ",F8.2," seconds")')timer_get_value(3)
  call timestamp(151)
  call flushifc(151)
endif

! generate G+k vectors for entire BZ (this is required to compute 
!   wave-functions at each k-point)
allocate(vgklnr(3,ngkmax,nkptnrloc))
allocate(vgkcnr(3,ngkmax,nkptnrloc))
allocate(gknr(ngkmax,nkptnrloc))
allocate(tpgknr(2,ngkmax,nkptnrloc))
allocate(ngknr(nkptnrloc))
allocate(sfacgknr(ngkmax,natmtot,nkptnrloc))
allocate(igkignr(ngkmax,nkptnrloc))
allocate(ylmgknr(lmmaxvr,ngkmax,nkptnrloc))
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  call gengpvec(vklnr(1,ik),vkcnr(1,ik),ngknr(ikloc),igkignr(1,ikloc), &
    vgklnr(1,1,ikloc),vgkcnr(1,1,ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc))
  call gensfacgp(ngknr(ikloc),vgkcnr(1,1,ikloc),ngkmax,sfacgknr(1,1,ikloc))
  do ig=1,ngknr(ikloc)
    call genylm(lmaxvr,tpgknr(1,ig,ikloc),ylmgknr(1,ig,ikloc))
  enddo
enddo
allocate(wfsvmtloc(lmmaxvr,nrfmax,natmtot,nspinor,nstsv,nkptnrloc))
allocate(wfsvitloc(ngkmax,nspinor,nstsv,nkptnrloc))
allocate(wfsvcgloc(ngkmax,nspinor,nstsv,nkptnrloc))  
allocate(evecfvloc(nmatmax,nstfv,nspnfv,nkptnrloc))
allocate(evecsvloc(nstsv,nstsv,nkptnrloc))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
if (wproc) then
  sz=lmmaxvr*nrfmax*natmtot*nstsv*nspinor
  sz=sz+ngkmax*nstsv*nspinor
  sz=sz+nmatmax*nstfv*nspnfv
  sz=sz+nstsv*nstsv
  sz=16*sz*nkptnrloc/1024/1024
  write(151,*)
  write(151,'("Size of wave-function arrays (MB) : ",I6)')sz
  write(151,*)
  write(151,'("Reading eigen-vectors")')
  call flushifc(151)
endif
call timer_start(1,reset=.true.)
! read eigen-vectors
if (mpi_grid_side(dims=(/dim_k/))) then
  do i=0,mpi_grid_size(dim_k)-1
    if (i.eq.mpi_grid_x(dim_k)) then
      do ikloc=1,nkptnrloc
        ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
        call getevecfv(vklnr(1,ik),vgklnr(1,1,ikloc),evecfvloc(1,1,1,ikloc))
        call getevecsv(vklnr(1,ik),evecsvloc(1,1,ikloc))
      enddo !ikloc
    endif
    if (.not.parallel_read) call mpi_grid_barrier(dims=(/dim_k/))
  enddo
endif !mpi_grid_side(dims=(/dim_k/)
call mpi_grid_barrier
call mpi_grid_bcast(evecfvloc(1,1,1,1),nmatmax*nstfv*nspnfv*nkptnrloc,&
  dims=(/dim2/))
call mpi_grid_bcast(evecsvloc(1,1,1),nstsv*nstsv*nkptnrloc,&
  dims=(/dim2/))
! transform eigen-vectors
wfsvmtloc=zzero
wfsvitloc=zzero
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
! get apw coeffs 
  call match(ngknr(ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc),        &
    sfacgknr(1,1,ikloc),apwalm)
! generate wave functions in muffin-tins
  call genwfsvmt(lmaxvr,lmmaxvr,ngknr(ikloc),evecfvloc(1,1,1,ikloc), &
    evecsvloc(1,1,ikloc),apwalm,wfsvmtloc(1,1,1,1,1,ikloc))
! generate wave functions in interstitial
  call genwfsvit(ngknr(ikloc),evecfvloc(1,1,1,ikloc), &
    evecsvloc(1,1,ikloc),wfsvitloc(1,1,1,ikloc))
enddo !ikloc
call timer_stop(1)
if (wproc) then
  write(151,'("Done in ",F8.2," seconds")')timer_get_value(1)
  call timestamp(151)
  call flushifc(151)
endif
! generate Wannier function expansion coefficients
if (wannier) then
  call timer_start(1,reset=.true.)
  if (allocated(wann_c)) deallocate(wann_c)
  allocate(wann_c(nwann,nstsv,nkptnrloc))
  if (allocated(wann_unkmt)) deallocate(wann_unkmt)
  allocate(wann_unkmt(lmmaxvr,nrfmax,natmtot,nspinor,nwann,nkptnrloc))
  if (allocated(wann_unkit)) deallocate(wann_unkit)
  allocate(wann_unkit(ngkmax,nspinor,nwann,nkptnrloc))
  wann_unkmt=zzero
  wann_unkit=zzero
  if (wproc) then
    write(151,*)
    write(151,'("Generating Wannier functions")')
    call flushifc(151)
  endif !wproc1
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    call genwann_c(ik,lr_evalsvnr(1,ik),wfsvmtloc(1,1,1,1,1,ikloc),&
      wann_c(1,1,ikloc))  
    do n=1,nwann
      do j=1,nstsv
        wann_unkmt(:,:,:,:,n,ikloc)=wann_unkmt(:,:,:,:,n,ikloc) + &
          wfsvmtloc(:,:,:,:,j,ikloc)*wann_c(n,j,ikloc)
        wann_unkit(:,:,n,ikloc)=wann_unkit(:,:,n,ikloc) + &
          wfsvitloc(:,:,j,ikloc)*wann_c(n,j,ikloc)
      enddo
    enddo
  enddo !ikloc
endif !wannier
! after optinal band disentanglement we can finally synchronize all eigen-values
!   and compute band occupation numbers 
call mpi_grid_reduce(lr_evalsvnr(1,1),nstsv*nkptnr,dims=(/dim_k/),all=.true.)
allocate(lr_occsvnr(nstsv,nkptnr))
call occupy2(nkptnr,wkptnr,lr_evalsvnr,lr_occsvnr)
deallocate(wfsvmtloc)
deallocate(wfsvitloc)
deallocate(wfsvcgloc)  
deallocate(evecfvloc)
deallocate(evecsvloc)
deallocate(apwalm)
deallocate(wann_c)


ntr=(2*trmax+1)**3
allocate(vtrl(3,ntr))
allocate(ivtit(-trmax:trmax,-trmax:trmax,-trmax:trmax))
n=0
do i1=-trmax,trmax
  do i2=-trmax,trmax
    do i3=-trmax,trmax
      n=n+1
      vtrl(:,n)=(/i1,i2,i3/)
      ivtit(i1,i2,i3)=n
    enddo
  enddo
enddo
ntrloc=mpi_grid_map(ntr,dim2)
if (wproc) then
  write(151,*)
  write(151,'("Number of translations : ",I4)')ntr
endif

allocate(wanmt(lmmaxvr,nrmtmax,natmtot,nspinor,ntrloc,nwann))
allocate(wanir(ngrtot,nspinor,ntrloc,nwann))
if (wproc) then
  sz=lmmaxvr*nrmtmax*natmtot+ngrtot
  sz=16*sz*nspinor*nwann*ntrloc/1024/1024
  write(151,*)
  write(151,'("Size of real-space Wannier functions arrays (MB) : ",I6)')sz
  write(151,*)
  call flushifc(151)
endif

call timer_reset(1)
call timer_reset(2)
allocate(wanmt0(lmmaxvr,nrmtmax,natmtot,nspinor,nwann))
allocate(wanir0(ngrtot,nspinor,nwann))
do itrloc=1,ntrloc
  itr=mpi_grid_map(ntr,dim2,loc=itrloc)
  call gen_wann_func(vtrl(1,itr),ngknr,vgkcnr,wanmt0,wanir0)
  do n=1,nwann
    wanmt(:,:,:,:,itrloc,n)=wanmt0(:,:,:,:,n)
    wanir(:,:,itrloc,n)=wanir0(:,:,n)
  enddo !n
enddo !itr
deallocate(wanmt0,wanir0)
if (wproc) then
  write(151,*)
  write(151,'("MT part : ",F8.3)')timer_get_value(1)
  write(151,'("IT part : ",F8.3)')timer_get_value(2)
endif

! calculate overlap matrix
allocate(ovlm(nwann,nwann))
if (mpi_grid_root(dims=(/dim_k/))) then
  ovlm=zzero
  do n1=1,nwann
    do n2=1,nwann
      ovlm(n1,n2)=inner_product(ntr,ntrloc,trmax,vtrl,ivtit,(/1,0,0/),&
        wanmt(1,1,1,1,1,n1),wanir(1,1,1,n1),wanmt(1,1,1,1,1,n2),&
        wanir(1,1,1,n2))
    enddo
  enddo
endif
if (wproc) then
  write(151,*)
  write(151,'("Overlap matrix")')
  write(151,'("  Real part : ")')  
  do n1=1,nwann
    write(151,'(2X,255F12.6)')(dreal(ovlm(n1,n2)),n2=1,nwann)
  enddo
  write(151,'("  Image part : ")')  
  do n1=1,nwann
    write(151,'(2X,255F12.6)')(dimag(ovlm(n1,n2)),n2=1,nwann)
  enddo
endif
if (wproc) close(151)
return
end

