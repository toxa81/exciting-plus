subroutine sic
use modmain
use mod_lfa
use mod_nrkp
implicit none

integer i1,i2,i3,n
integer ik,ikloc,ik1,j,ig,sz,i,isym
integer n1,n2,ispn
integer itr,ntrloc,itrloc

complex(8), allocatable :: wanmt(:,:,:,:,:,:)
complex(8), allocatable :: wanir(:,:,:,:)
complex(8), allocatable :: wanmt0(:,:,:,:,:)
complex(8), allocatable :: wanir0(:,:,:)

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

call lfa_init(0)
call genwfnr(151)
! deallocate unnecessary wave-functions
deallocate(wfsvmtloc)
deallocate(wfsvitloc)
deallocate(evecfvloc)
deallocate(evecsvloc)
deallocate(wann_c)

ntrloc=mpi_grid_map(ntr,dim2)
if (wproc) then
  write(151,*)
  write(151,'("Number of translations : ",I4)')ntr
endif

allocate(wanmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor,nwann))
allocate(wanir(ngrtot,ntrloc,nspinor,nwann))
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
  call gen_wann_func(vtl(1,itr),ngknr,vgkcnr,wanmt0,wanir0)
  do ispn=1,nspinor
    do n=1,nwann
      wanmt(:,:,:,itrloc,ispn,n)=wanmt0(:,:,:,ispn,n)
      wanir(:,itrloc,ispn,n)=wanir0(:,ispn,n)
    enddo !n
  enddo !ispn
enddo !itr
deallocate(wanmt0,wanir0)
if (wproc) then
  write(151,*)
  write(151,'("MT part : ",F8.3)')timer_get_value(1)
  write(151,'("IT part : ",F8.3)')timer_get_value(2)
endif

! calculate overlap matrix
allocate(ovlm(nwann,nwann))
ovlm=zzero
if (mpi_grid_root(dims=(/dim_k/))) then
  ovlm=zzero
  do n1=1,nwann
    do n2=1,nwann
      do ispn=1,nspinor
        ovlm(n1,n2)=ovlm(n1,n2)+lfa_dotp(.true.,(/0,0,0/),wanmt(1,1,1,1,ispn,n1),&
          wanir(1,1,ispn,n1),wanmt(1,1,1,1,ispn,n2),wanir(1,1,ispn,n2))
      enddo
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

