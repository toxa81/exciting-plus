#ifdef _HDF5_
subroutine response
use modmain
use hdf5
implicit none

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
complex(8), allocatable :: pmat(:,:,:,:)

integer i,j,n,ik,ikloc,ik1,isym,idx0,ivq1,ivq2,iq,ig
integer sz,nvq0loc,i1,i2,i3
character*100 qnm
real(8) w2
logical lgamma,lpmat,wproc1
logical, external :: wann_diel

! comment: after all new implementations (response in WF, cRPA,
!   band disentanglement) the code has become ugly and unreadable;
!   it should be refactored; new hdf5 and mpi_grid interfaces are
!   "a must"
!
! after some refactoring the code is still ugly; do something with wannir
!   response
!
! typical execution patterns
!  I) compute and save ME (task 400), read ME, compute and save chi0 (task 401),
!     read chi0 and compute chi (task 402)
!  II) the same + Wannier channels decomposition 
!  III) same as I) and II) but without saving matrix elements
!  IV) constraind RPA: compute ME, compute chi0 -> chi -> screened W and U 
!
! New task list:
!   400 - compute and write ME
!   401 - compute and write chi0
!   402 - compute and write chi
!   403 - compute ME, compute and write chi0, compute and write chi OR
!         compute ME, compute chi0 chi and screened U, sum over q; behaviour is
!         contorlled by crpa
!
! MPI grid for tasks:
!   400 (matrix elements) : (1) k-points x (2) G-vectors or interband 
!                                           transitions x (3) q-points 
!   401 (chi0) : (1) k-points x (2) interband transition x (3) q-points 
!   402 (chi) : (1) energy mesh x (2) number of fxc kernels x (3) q-points
!
! todo: more comments!!!

if (lrtype.eq.1.and..not.spinpol) then
  write(*,*)
  write(*,'("Error(response): can''t do magnetic response for unpolarized ground state")')
  write(*,*)
  call pstop
endif
if (nvq0.eq.0) then
  write(*,*)
  write(*,'("Error(response): no q-vectors")')
  write(*,*)
  call pstop
endif
if (crpa.and.task.ne.403) then
  write(*,*)
  write(*,'("Error(response): cRPA must be run with task=403")')
  write(*,*)
  call pstop
endif

if (.not.wannier) then
  wannier_chi0_chi=.false.
  lwannopt=.false.
  crpa=.false.
endif
lpmat=.false.
!if (lwannopt.or.crpa) lpmat=.true.
do j=1,nvq0
  if (ivq0m_list(1,j).eq.0.and.ivq0m_list(2,j).eq.0.and.ivq0m_list(3,j).eq.0) lpmat=.true.
enddo
! set the switch to write matrix elements
write_megq_file=.true.
if (task.eq.403) write_megq_file=.false.
write_chi0_file=.true.
if (crpa) write_chi0_file=.false.

! set the switch to compute screened W matrix in task 402
screened_w=.false.
!if (crpa) screened_w=.true.

wannier_megq=.false.
if (crpa.or.wannier_chi0_chi) wannier_megq=.true.

! this is enough for matrix elements
lmaxvr=5

! initialise universal variables
call init0
call init1

if (.not.mpi_grid_in()) return

! for constrained RPA all q-vectors in BZ are required 
lgamma=.false.
if (crpa) then
  if (allocated(ivq0m_list)) deallocate(ivq0m_list)
  if (lgamma) then
    nvq0=nkptnr
  else
    nvq0=nkptnr-1 
  endif
  allocate(ivq0m_list(3,nvq0))
  j=0
  do i1=0,ngridk(1)-1
    do i2=0,ngridk(2)-1
      do i3=0,ngridk(3)-1
        if (.not.(i1.eq.0.and.i2.eq.0.and.i3.eq.0.and..not.lgamma)) then
          j=j+1
          ivq0m_list(:,j)=(/i1,i2,i3/)
        endif
      enddo
    enddo
  enddo
endif
! todo: put warnings to output
if (crpa) then
  maxomega=0.d0
  domega=1.d0
  nfxca=1
  fxca0=0.d0
  fxca1=0.d0
endif

! necessary calls before generating Bloch wave-functions 
if (task.eq.400.or.task.eq.403) then
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
endif
! create q-directories
if (mpi_grid_root()) then
  do iq=1,nvq0
    call qname(ivq0m_list(:,iq),qnm)
    call system("mkdir -p "//trim(qnm))
  enddo
endif
call mpi_grid_barrier()

wproc1=.false.
if (mpi_grid_root()) then
  wproc1=.true.
  if (task.eq.400) open(151,file='RESPONSE_ME.OUT',form='formatted',status='replace')
  if (task.eq.401) open(151,file='RESPONSE_CHI0.OUT',form='formatted',status='replace')
  if (task.eq.402) open(151,file='RESPONSE_CHI.OUT',form='formatted',status='replace')
  if (task.eq.403) open(151,file='RESPONSE.OUT',form='formatted',status='replace')  
  call timestamp(151)
endif
if (wproc1) then
  write(151,'("Total number of q-vectors        : ",I6)')nvq0
  write(151,'("Total number of processors       : ",I6)')nproc
  write(151,'("MPI grid size                    : ",3I6)')mpi_grid_size
  if (nproc.gt.1) then
    write(151,'("Parallel file reading            : ",L1)')parallel_read
  endif
  if (nproc.gt.1) then
    write(151,'("Parallel file writing            : ",L1)')parallel_write
  endif
  write(151,'("Wannier functions                : ",L1)')wannier
  write(151,'("Response in Wannier basis        : ",L1)')wannier_chi0_chi
  write(151,'("Constrained RPA                  : ",L1)')crpa
  write(151,'("Matrix elements in Wannier basis : ",L1)')wannier_megq
  write(151,'("Write matrix elements            : ",L1)')write_megq_file
  write(151,'("Write chi0                       : ",L1)')write_chi0_file  
  call flushifc(151)
endif

if (task.eq.400.or.task.eq.403) then
! get energies of states in reduced part of BZ
  call timer_start(3,reset=.true.)
  if (wproc1) then
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
  if (wproc1) then
    write(151,'("Done in ",F8.2," seconds")')timer_get_value(3)
    call timestamp(151)
    call flushifc(151)
  endif
endif

! generate wave-functions
if (task.eq.400.or.task.eq.403) then
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
  if (lpmat) then
    allocate(pmat(3,nstsv,nstsv,nkptnrloc))
  endif
  if (wproc1) then
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
    dims=(/dim2,dim3/))
  call mpi_grid_bcast(evecsvloc(1,1,1),nstsv*nstsv*nkptnrloc,&
    dims=(/dim2,dim3/))
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
! generate <e^{-i(G+k)x}|\psi>
!    call genwfsvcg(ngknr(ikloc),igkignr(1,ikloc),gknr(1,ikloc), &
!      ylmgknr(1,1,ikloc),sfacgknr(1,1,ikloc),wfsvmtloc(1,1,1,1,1,ikloc), &
!        wfsvitloc(1,1,1,ikloc),wfsvcgloc(1,1,1,ikloc))
!        call pstop
    if (lpmat) then
      call genpmat(ngknr(ikloc),igkignr(1,ikloc),vgkcnr(1,1,ikloc),&
        apwalm,evecfvloc(1,1,1,ikloc),evecsvloc(1,1,ikloc),pmat(1,1,1,ikloc))
    endif
  enddo !ikloc
  call timer_stop(1)
  if (wproc1) then
    write(151,'("Done in ",F8.2," seconds")')timer_get_value(1)
    call timestamp(151)
    call flushifc(151)
  endif
! generate Wannier function expansion coefficients
  if (wannier_megq) then
    call timer_start(1,reset=.true.)
    if (allocated(wann_c)) deallocate(wann_c)
! use first nkptnrloc points to store wann_c(k) and second nkptnrloc points
!   to store wann_c(k+q)
    allocate(wann_c(nwann,nstsv,2*nkptnrloc))
    if (wproc1) then
      write(151,*)
      write(151,'("Generating Wannier functions")')
      call flushifc(151)
    endif !wproc1
    do ikloc=1,nkptnrloc
      ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
      call genwann_c(ik,lr_evalsvnr(1,ik),wfsvmtloc(1,1,1,1,1,ikloc),&
        wann_c(1,1,ikloc))  
      if (ldisentangle) then
! disentangle bands
        call disentangle(lr_evalsvnr(1,ik),wann_c(1,1,ikloc),evecsvloc(1,1,ikloc))
! recompute wave functions
! get apw coeffs 
        call match(ngknr(ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc),        &
          sfacgknr(1,1,ikloc),apwalm)
! generate wave functions in muffin-tins
        call genwfsvmt(lmaxvr,lmmaxvr,ngknr(ikloc),evecfvloc(1,1,1,ikloc), &
          evecsvloc(1,1,ikloc),apwalm,wfsvmtloc(1,1,1,1,1,ikloc))
! generate wave functions in interstitial
        call genwfsvit(ngknr(ikloc),evecfvloc(1,1,1,ikloc), &
          evecsvloc(1,1,ikloc),wfsvitloc(1,1,1,ikloc))       
      endif
    enddo !ikloc
  endif !wannier
! after optinal band disentanglement we can finally synchronize all eigen-values
!   and compute band occupation numbers 
  call mpi_grid_reduce(lr_evalsvnr(1,1),nstsv*nkptnr,dims=(/dim_k/),all=.true.)
  allocate(lr_occsvnr(nstsv,nkptnr))
  call occupy2(nkptnr,wkptnr,lr_evalsvnr,lr_occsvnr)
  if (wannier_megq) then
! calculate Wannier function occupancies 
    wann_occ=0.d0
    do n=1,nwann
      do ikloc=1,nkptnrloc
        ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
        do j=1,nstsv
          w2=dreal(dconjg(wann_c(n,j,ikloc))*wann_c(n,j,ikloc))
          wann_occ(n)=wann_occ(n)+w2*lr_occsvnr(j,ik)/nkptnr
        enddo
      enddo
    enddo
    call mpi_grid_reduce(wann_occ(1),nwann,dims=(/dim_k/),all=.true.)
    if (wproc1) then
      write(151,'("  Wannier function occupation numbers : ")')
      do n=1,nwann
        write(151,'("    n : ",I4,"  occ : ",F8.6)')n,wann_occ(n)
      enddo
    endif
    if (wproc1) then
      write(151,'("  Dielectric Wannier functions : ",L1)')wann_diel()
    endif
    call timer_stop(1)
    if (wproc1) then
      write(151,'("Done in ",F8.2," seconds")')timer_get_value(1)
      call timestamp(151)
      call flushifc(151)
    endif
  endif !wannier
  deallocate(apwalm)
  deallocate(vgklnr)
  deallocate(vgkcnr)
  deallocate(gknr)
  deallocate(tpgknr)
  deallocate(sfacgknr)
endif !task.eq.400.or.task.eq.403

! distribute q-vectors along 3-rd dimention
nvq0loc=mpi_grid_map(nvq0,dim3,offs=idx0)
ivq1=idx0+1
ivq2=idx0+nvq0loc

!-----------------------------------------!
!    task 400: compute matrix elements    !
!-----------------------------------------!
if (task.eq.400) then
! calculate matrix elements
  call timer_start(10,reset=.true.)
  if (wproc1) call timestamp(151,txt='start genmegq')
  do iq=ivq1,ivq2
    call genmegq(ivq0m_list(1,iq),wfsvmtloc,wfsvitloc,ngknr, &
      igkignr,pmat)
  enddo
  call timer_stop(10)
  if (wproc1) call timestamp(151,txt='stop genmegq')
  if (wproc1) then
    write(151,*)
    write(151,'("Total time for matrix elements : ",F8.2," seconds")')timer_get_value(10)
    call flushifc(151)
  endif
endif

!------------------------------!
!    task 401: compute chi0    !
!------------------------------!
if (task.eq.401) then
! calculate chi0
  call timer_start(11,reset=.true.)
  if (wproc1) call timestamp(151,txt='start genchi0')
  do iq=ivq1,ivq2
    call genchi0(ivq0m_list(1,iq))
  enddo
  call timer_stop(11)
  if (wproc1) call timestamp(151,txt='stop genchi0')
  if (wproc1) then
    write(151,*)
    write(151,'("Total time for chi0 : ",F8.2," seconds")')timer_get_value(11)
    call flushifc(151)
  endif
endif

!-----------------------------!
!    task 402: compute chi    !
!-----------------------------!
if (task.eq.402) then
! calculate chi
  call timer_start(12,reset=.true.)
  do iq=ivq1,ivq2
    call genchi(ivq0m_list(1,iq))
  enddo
  call timer_stop(12)
  if (wproc1) then
    write(151,*)    
    write(151,'("Total time for chi : ",F8.2," seconds")')timer_get_value(12)
    call flushifc(151)
  endif
endif

!------------------------------------------!
!    task 403: compute me, chi0 and chi    !
!------------------------------------------!
if (task.eq.403) then
  do iq=ivq1,ivq2
    call genmegq(ivq0m_list(1,iq),wfsvmtloc,wfsvitloc,ngknr, &
      igkignr,pmat)
    call genchi0(ivq0m_list(1,iq))
    if (.not.crpa) call genchi(ivq0m_list(1,iq))
  enddo
  call mpi_grid_barrier()
  if (crpa.and.mpi_grid_root()) call response_u
endif

if (wproc1) close(151)

if (task.eq.400.and.lpmat) deallocate(pmat)

if (task.eq.400.or.task.eq.403) then
  deallocate(wfsvmtloc)
  deallocate(wfsvitloc)
  deallocate(evecfvloc)
  deallocate(evecsvloc)
  deallocate(ngknr)
  deallocate(igkignr)
  deallocate(lr_occsvnr)
  deallocate(lr_evalsvnr)   
  if (wannier_megq) deallocate(wann_c)
endif

return
end
#endif