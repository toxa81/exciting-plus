#ifdef _HDF5_
subroutine response
use modmain
#ifdef _MPI_
use mpi
#endif
use hdf5
implicit none

integer, allocatable :: ngknr(:)
integer, allocatable :: igkignr(:,:)
real(8), allocatable :: vgklnr(:,:,:)
real(8), allocatable :: vgkcnr(:,:,:)
real(8), allocatable :: gknr(:,:)
real(8), allocatable :: tpgknr(:,:,:)
complex(8), allocatable :: sfacgknr(:,:,:)

real(8), allocatable :: occsvnr(:,:)
real(8), allocatable :: evalsvnr(:,:)
complex(8), allocatable :: wfsvitloc(:,:,:,:)
complex(8), allocatable :: wfsvmtloc(:,:,:,:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: pmat(:,:,:,:)

integer i,j,n,ik,ikloc,istsv,ik1,isym,idx0,bs,ivq1,ivq2,iq
integer sz,iint,iw
integer i1,i2,i3
character*100 fname,qnm
character*3 c3
real(8) w2
logical lgamma,lpmat
integer, external :: iknrglob2
logical, external :: root_cart
logical, external :: in_cart
logical, external :: wann_diel

! comment: after all new implementations (response in WF, cRPA,
!   band disentanglement) the code has become ugly and unreadable;
!   it should be refactored; new hdf5 and mpi_grid interfaces are
!   "a must"

if (lrtype.eq.1.and..not.spinpol) then
  write(*,*)
  write(*,'("Error(response): can''t do magnetic response for unpolarized ground state")')
  write(*,*)
  call pstop
endif
if (.not.wannier) then
  lwannresp=.false.
  lwannopt=.false.
  crpa=.false.
endif
lpmat=.false.
!if (lwannopt.or.crpa) lpmat=.true.
do j=1,nvq0
  if (ivq0m_list(1,j).eq.0.and.ivq0m_list(2,j).eq.0.and.ivq0m_list(3,j).eq.0) lpmat=.true.
enddo
! this is enough for matrix elements
lmaxvr=4
! initialise universal variables
call init0
call init1
! for constrained RPA all q-vectors in BZ are required 
lgamma=.false.
if (crpa.and..false.) then
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
if (crpa) then
  maxomega=0.d0
  domega=1.d0
endif
! parallel init; MPI grid for tasks:
!   400 (matrix elements) : (1) k-points x (2) G-vectors x (3) q-points 
!   401 (chi0) : (1) k-points x (2) interband transition x (3) q-points 
!   402 (chi) : (1) energy mesh x (2) number of fxc kernels x (3) q-points
call pinit_cart
! necessary calls before generating Bloch wave-functions 
if (task.eq.400) then
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
  if (iproc.eq.0) call readfermi
  call dsync(efermi,1,.false.,.true.)
endif
! create q-directories
if (iproc.eq.0) then
  do iq=1,nvq0
    call qname(ivq0m_list(:,iq),qnm)
    call system("mkdir -p "//trim(qnm))
  enddo
endif
call barrier(comm_world)
! work only with processors in the grid
if (in_cart()) then
  wproc=.false.
  if (root_cart((/1,1,1/))) then
    wproc=.true.
    if (task.eq.400) open(151,file='RESPONSE_ME.OUT',form='formatted',status='replace')
    if (task.eq.401) open(151,file='RESPONSE_CHI0.OUT',form='formatted',status='replace')
    if (task.eq.402) open(151,file='RESPONSE_CHI.OUT',form='formatted',status='replace')
    if (task.eq.403) open(151,file='RESPONSE_U.OUT',form='formatted',status='replace')
  endif
  if (wproc) then
    write(151,'("Running on ",I8," proc.")')nproc
#ifdef _PIO_
    if (nproc.gt.1) then
      write(151,'("Using parallel I/O")')
    endif
#endif
    write(151,'("MPI grid size : ",3I6)')mpi_dims
    write(151,'("Wannier functions switch : ",L1)')wannier
    write(151,'("Response in local basis  : ",L1)')lwannresp
    call flushifc(151)
  endif
  
! distribute k-points over 1-st dimension of the grid
  if (allocated(nkptnrloc)) deallocate(nkptnrloc)
  allocate(nkptnrloc(0:mpi_dims(1)-1))
  if (allocated(ikptnrloc)) deallocate(ikptnrloc)
  allocate(ikptnrloc(0:mpi_dims(1)-1,2))
  call splitk(nkptnr,mpi_dims(1),nkptnrloc,ikptnrloc)
  nkptnr_loc=nkptnrloc(mpi_x(1))
  
  if (task.eq.400) then
! get energies of states in reduced part of BZ
!    allocate(occsvnr(nstsv,nkptnr))
    call timer_start(3)
    if (wproc) then
      write(151,*)
      write(151,'("Reading energies of states")')
      call flushifc(151)
    endif
    if (root_cart((/1,1,1/))) then
! read from IBZ
      do ik=1,nkpt
        call getevalsv(vkl(1,ik),evalsv(1,ik))
      enddo
    endif
    call d_bcast_cart(comm_cart,evalsv,nstsv*nkpt)
    allocate(evalsvnr(nstsv,nkptnr))
    evalsvnr=0.d0
    do ikloc=1,nkptnr_loc
      ik=iknrglob2(ikloc,mpi_x(1))
      call findkpt(vklnr(1,ik),isym,ik1) 
      evalsvnr(:,ik)=evalsv(:,ik1)
    enddo
    call timer_stop(3)
    if (wproc) then
      write(151,'("Done in ",F8.2," seconds")')timer(3,2)
      call flushifc(151)
    endif
  endif
  
  if (task.eq.400) then
! generate G+k vectors for entire BZ (this is required to compute 
!   wave-functions at each k-point)
    allocate(vgklnr(3,ngkmax,nkptnr_loc))
    allocate(vgkcnr(3,ngkmax,nkptnr_loc))
    allocate(gknr(ngkmax,nkptnr_loc))
    allocate(tpgknr(2,ngkmax,nkptnr_loc))
    allocate(ngknr(nkptnr_loc))
    allocate(sfacgknr(ngkmax,natmtot,nkptnr_loc))
    allocate(igkignr(ngkmax,nkptnr_loc))
    do ikloc=1,nkptnr_loc
      ik=iknrglob2(ikloc,mpi_x(1))
      call gengpvec(vklnr(1,ik),vkcnr(1,ik),ngknr(ikloc),igkignr(1,ikloc), &
        vgklnr(1,1,ikloc),vgkcnr(1,1,ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc))
      call gensfacgp(ngknr(ikloc),vgkcnr(1,1,ikloc),ngkmax,sfacgknr(1,1,ikloc))
    enddo
    allocate(wfsvmtloc(lmmaxvr,nrfmax,natmtot,nspinor,nstsv,nkptnr_loc))
    allocate(wfsvitloc(ngkmax,nspinor,nstsv,nkptnr_loc))
    allocate(evecfvloc(nmatmax,nstfv,nspnfv,nkptnr_loc))
    allocate(evecsvloc(nstsv,nstsv,nkptnr_loc))
    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
    if (lpmat) then
      allocate(pmat(3,nstsv,nstsv,nkptnr_loc))
    endif
    if (wproc) then
      sz=lmmaxvr*nrfmax*natmtot*nstsv*nspinor
      sz=sz+ngkmax*nstsv*nspinor
      sz=sz+nmatmax*nstfv*nspnfv
      sz=sz+nstsv*nstsv
      sz=16*sz*nkptnr_loc/1024/1024
      write(151,*)
      write(151,'("Size of wave-function arrays (MB) : ",I6)')sz
      write(151,*)
      write(151,'("Reading eigen-vectors")')
      call flushifc(151)
    endif
    call timer_start(1)
! read and transform eigen-vectors
    wfsvmtloc=zzero
    wfsvitloc=zzero
    if (root_cart((/0,1,1/))) then
#ifndef _PIO_
      do i=0,mpi_dims(1)-1
        if (i.eq.mpi_x(1)) then
#endif
          do ikloc=1,nkptnr_loc
            ik=iknrglob2(ikloc,mpi_x(1))
            call getevecfv(vklnr(1,ik),vgklnr(1,1,ikloc),evecfvloc(1,1,1,ikloc))
            call getevecsv(vklnr(1,ik),evecsvloc(1,1,ikloc))
! get apw coeffs 
            call match(ngknr(ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc),        &
              sfacgknr(1,1,ikloc),apwalm)
! generate wave functions in muffin-tins
            call genwfsvmt(lmaxvr,lmmaxvr,ngknr(ikloc),evecfvloc(1,1,1,ikloc), &
              evecsvloc(1,1,ikloc),apwalm,wfsvmtloc(1,1,1,1,1,ikloc))
! generate wave functions in interstitial
            call genwfsvit(ngknr(ikloc),evecfvloc(1,1,1,ikloc), &
              evecsvloc(1,1,ikloc),wfsvitloc(1,1,1,ikloc))
            if (lpmat) then
              call genpmat(ngknr(ikloc),igkignr(1,ikloc),vgkcnr(1,1,ikloc),&
                apwalm,evecfvloc(1,1,1,ikloc),evecsvloc(1,1,ikloc),pmat(1,1,1,ikloc))
            endif
          enddo !ikloc
#ifndef _PIO_
        endif
        call barrier(comm_cart_100)
      enddo
#endif
    endif !root_cart
    call barrier(comm_cart_110)
    call d_bcast_cart(comm_cart_011,wfsvmtloc, &
      lmmaxvr*nrfmax*natmtot*nspinor*nstsv*nkptnr_loc*2)
    call d_bcast_cart(comm_cart_011,wfsvitloc, &
      ngkmax*nspinor*nstsv*nkptnr_loc*2)
    call d_bcast_cart(comm_cart_011,evecfvloc, &
      nmatmax*nstfv*nspnfv*nkptnr_loc*2)
    call d_bcast_cart(comm_cart_011,evecsvloc, &
      nstsv*nstsv*nkptnr_loc*2)
    call timer_stop(1)
    if (wproc) then
      write(151,'("Done in ",F8.2," seconds")')timer(1,2)
      call flushifc(151)
    endif
    call timer_reset(1)
    call timer_start(1)
! generate Wannier function expansion coefficients
    if (wannier) then
      if (allocated(wann_c)) deallocate(wann_c)
      allocate(wann_c(nwann,nstsv,nkptnr_loc))
      if (wproc) then
        write(151,*)
        write(151,'("Generating Wannier functions")')
        call flushifc(151)
      endif !wproc
      do ikloc=1,nkptnr_loc
        ik=iknrglob2(ikloc,mpi_x(1))
        call genwann_c(ik,evalsvnr(1,ik),wfsvmtloc(1,1,1,1,1,ikloc),wann_c(1,1,ikloc))
        if (ldisentangle) then
! disentangle bands
          call disentangle(evalsvnr(1,ik),wann_c(1,1,ikloc),evecsvloc(1,1,ikloc))
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
    call timer_stop(1)
    if (wproc) then
      write(151,'("Done in ",F8.2," seconds")')timer(1,2)
      call flushifc(151)
    endif
! after optinal band disentanglement we can finally synchronize all eigen-values
!   and compute band occupation numbers 
    call d_reduce_cart(comm_cart_100,.true.,evalsvnr,nstsv*nkptnr)
    allocate(occsvnr(nstsv,nkptnr))
    call occupy2(nkptnr,wkptnr,evalsvnr,occsvnr)
    if (wannier) then
! calculate Wannier function occupancies 
      wann_occ=0.d0
      do n=1,nwann
        do ikloc=1,nkptnr_loc
          ik=iknrglob2(ikloc,mpi_x(1))
          do j=1,nstsv
            w2=dreal(dconjg(wann_c(n,j,ikloc))*wann_c(n,j,ikloc))
            wann_occ(n)=wann_occ(n)+w2*occsvnr(j,ik)/nkptnr
          enddo
        enddo
      enddo
      call d_reduce_cart(comm_cart_100,.true.,wann_occ,nwann)
      if (wproc) then
        write(151,'("  Wannier function occupation numbers : ")')
        do n=1,nwann
          write(151,'("    n : ",I4,"  occ : ",F8.6)')n,wann_occ(n)
        enddo
      endif
      if (wproc) then
        write(151,'("  Dielectric Wannier functions : ",L1)')wann_diel()
      endif
    endif !wannier
    call timer_stop(1)
    if (wproc) then
      write(151,'("Done in ",F8.2," seconds")')timer(1,2)
      call flushifc(151)
    endif
    deallocate(apwalm)
    deallocate(vgklnr)
    deallocate(vgkcnr)
    deallocate(gknr)
    deallocate(tpgknr)
    deallocate(sfacgknr)
  endif !task.eq.400
  
! distribute q-vectors along 3-rd dimention 
  call idxbos(nvq0,mpi_dims(3),mpi_x(3)+1,idx0,bs)
  ivq1=idx0+1
  ivq2=idx0+bs
  
  if (task.eq.400) then
! calculate matrix elements
    call timer_start(10)
    do iq=ivq1,ivq2
      call response_me(ivq0m_list(1,iq),wfsvmtloc,wfsvitloc,ngknr, &
        igkignr,occsvnr,evalsvnr,pmat)
    enddo
    call timer_stop(10)
    if (wproc) then
      write(151,*)
      write(151,'("Total time for matrix elements : ",F8.2," seconds")')timer(10,2)
      call flushifc(151)
    endif
  endif
  
  if (task.eq.401) then
! calculate chi0
    call timer_start(11)
    do iq=ivq1,ivq2
      call response_chi0(ivq0m_list(1,iq))
    enddo
    call timer_stop(11)
    if (wproc) then
      write(151,*)
      write(151,'("Total time for chi0 : ",F8.2," seconds")')timer(11,2)
      call flushifc(151)
    endif
  endif
  
  if (task.eq.402) then
! calculate chi
    call timer_start(12)
    do iq=ivq1,ivq2
      call response_chi(ivq0m_list(1,iq))
    enddo
    call timer_stop(12)
    if (wproc) then
      write(151,*)    
      write(151,'("Total time for chi : ",F8.2," seconds")')timer(12,2)
      call flushifc(151)
    endif
  endif
  
  if (task.eq.403.and.crpa.and.root_cart((/1,1,1/))) then
    call response_u
  endif
    
  if (wproc) close(151)
  
  if (task.eq.400.and.lpmat) deallocate(pmat)
  
  if (task.eq.400) then
    deallocate(wfsvmtloc)
    deallocate(wfsvitloc)
    deallocate(evecfvloc)
    deallocate(evecsvloc)
    deallocate(ngknr)
    deallocate(igkignr)
    deallocate(occsvnr)
    deallocate(evalsvnr)   
    if (wannier) deallocate(wann_c)
  endif
endif !in_cart()

return
end
#endif