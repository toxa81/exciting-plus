#ifdef _HDF5_
subroutine response
use modmain
use hdf5
use mod_nrkp
implicit none
#ifdef _PAPI_
integer ierr
include 'f90papi.h'
real real_time,cpu_time,mflops
integer*8 fp_ins
#endif

integer i,j,n,ik,ikloc,ik1,isym,idx0,iq,ig
integer n1,nwloc
logical l1
integer sz,nvq0loc,iqloc,i1,i2,i3,ias,jas
character*100 qnm
real(8) w2,t1
logical lgamma,wproc1,lpmat
logical, external :: wann_diel

! Task list:
!   400 - response functions
!   401 - cRPA
!   402 - bare U
!
! MPI grid:
!   (matrix elements) : (1) k-points x (2) interband transitions x 
!                     x (3) q-points 
!   (response) : (1) energy mesh x (2) number of fxc kernels x (3) q-points
!

if (lrtype.eq.1.and..not.spinpol) then
  write(*,*)
  write(*,'("Error(response): can''t do magnetic response for unpolarized ground state")')
  write(*,*)
  call pstop
endif
if (nvq0.eq.0.and..not.(task.eq.401.or.task.eq.402)) then
  write(*,*)
  write(*,'("Error(response): no q-vectors")')
  write(*,*)
  call pstop
endif   
if (.not.wannier) wannier_chi0_chi=.false.
if (.not.spinpol) megqwan_afm=.false.
if (wannier_chi0_chi.or.task.eq.401.or.task.eq.402) wannier_megq=.true.

! this is enough for matrix elements
!lmaxvr=5

if (iproc.eq.0) call timestamp(6,'before init0')
! initialise universal variables
call init0
call init1
call mpi_world_barrier
if (iproc.eq.0) call timestamp(6,'after init1')
if (.not.mpi_grid_in()) return

! check if all q-vectors in BZ are required 
! TODO: put to a separate file
lgamma=.true.
if (task.eq.401.or.task.eq.402) then
  if (allocated(ivq0m_list)) deallocate(ivq0m_list)
  nvq0=nkptnr-1
  j=0  
  if (lgamma) then
    nvq0=nvq0+8
    j=8
  endif
  allocate(ivq0m_list(3,nvq0))
  ivq0m_list=0
  do i1=0,ngridk(1)-1
    do i2=0,ngridk(2)-1
      do i3=0,ngridk(3)-1
        if (.not.(i1.eq.0.and.i2.eq.0.and.i3.eq.0)) then
          j=j+1
          ivq0m_list(:,j)=(/i1,i2,i3/)
        endif
      enddo
    enddo
  enddo
  nfxca=1
  fxca0=0.d0
  fxca1=0.d0
  if (lgamma) call init_q0gamma
endif
! check if momentum matrix is required
lpmat=.false.
do j=1,nvq0
  if (ivq0m_list(1,j).eq.0.and.ivq0m_list(2,j).eq.0.and. &
      ivq0m_list(3,j).eq.0) then
    lpmat=.true.
  endif
enddo
! create q-directories
if (mpi_grid_root()) then
  call system("mkdir -p q")
  do iq=1,nvq0
    call getqdir(iq,ivq0m_list(:,iq),qnm)
    call system("mkdir -p "//trim(qnm))
  enddo
endif

wproc1=.false.
if (mpi_grid_root()) then
  wproc1=.true.
  open(151,file="RESPONSE.OUT",form="FORMATTED",status="REPLACE")
  call timestamp(151)
endif
wproc=wproc1
if (wproc1) then
  write(151,'("Total number of q-vectors        : ",I6)')nvq0
  write(151,'("Total number of processors       : ",I6)')nproc
  write(151,'("MPI grid size                    : ",3I6)')mpi_grid_size
  if (nproc.gt.1) then
    write(151,'("Parallel file reading            : ",L1)')parallel_read
  endif
  write(151,'("Wannier functions                : ",L1)')wannier
  write(151,'("Matrix elements in Wannier basis : ",L1)')wannier_megq
  write(151,'("Response in Wannier basis        : ",L1)')wannier_chi0_chi
  call flushifc(151)
endif

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
! generate wave-functions for entire BZ
call genwfnr(151,lpmat)
if (spinpol) then
  if (allocated(spinor_ud)) deallocate(spinor_ud)
  allocate(spinor_ud(2,nstsv,nkptnr))
  spinor_ud=0
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do j=1,nstsv
      t1=sum(abs(evecsvloc(1:nstfv,j,ikloc)))
      if (t1.gt.1d-10) spinor_ud(1,j,ik)=1
      t1=sum(abs(evecsvloc(nstfv+1:nstsv,j,ikloc)))
      if (t1.gt.1d-10) spinor_ud(2,j,ik)=1
    enddo
  enddo
  call mpi_grid_reduce(spinor_ud(1,1,1),2*nstsv*nkptnr,dims=(/dim_k/),all=.true.)
endif  

! TODO: put to separate file
if (wannier_megq) then
  call getnghbr(megqwan_maxdist)
  nmegqwanmax=0
  do n=1,nwann
    ias=iwann(1,n)
    do i=1,nnghbr(ias)
      do n1=1,nwann
        jas=iwann(1,n1)
        if (jas.eq.inghbr(1,i,ias)) then
          nmegqwanmax=nmegqwanmax+nwannias(jas)
        endif
      enddo
    enddo
  enddo
  allocate(imegqwan(5,nmegqwanmax))
  imegqwan=0
  nmegqwan=0   
  do n=1,nwann
    ias=iwann(1,n)
    do i=1,nnghbr(ias)
      do n1=1,nwann
        jas=iwann(1,n1)
        if (jas.eq.inghbr(1,i,ias)) then
          l1=.false.
! for integer occupancy numbers take only transitions between occupied and empty bands
          if (wann_diel().and.(abs(wann_occ(n)-wann_occ(n1)).gt.1d-8)) l1=.true.
! for fractional occupancies or cRPA calculation take all transitions
          if (.not.wann_diel().or.task.eq.401.or.task.eq.402) l1=.true.
          if (l1) then
            nmegqwan=nmegqwan+1
            imegqwan(1,nmegqwan)=n
            imegqwan(2,nmegqwan)=n1
            imegqwan(3:5,nmegqwan)=inghbr(3:5,i,ias)
          endif 
        endif
      enddo
    enddo
  enddo
  if (wproc1) then
    write(151,*)
    write(151,'("Number of Wannier transitions : ",I6)')nmegqwan
    write(151,'("List of Wannier transitions (n n1 T) ")')
    do i=1,nmegqwan
      write(151,'(I4,4X,I4,4X,3I3)')imegqwan(:,i)
    enddo
    call timestamp(151)
    call flushifc(151)
  endif    
  megqwan_tlim(1,1)=minval(imegqwan(3,:))
  megqwan_tlim(2,1)=maxval(imegqwan(3,:))
  megqwan_tlim(1,2)=minval(imegqwan(4,:))
  megqwan_tlim(2,2)=maxval(imegqwan(4,:))
  megqwan_tlim(1,3)=minval(imegqwan(5,:))
  megqwan_tlim(2,3)=maxval(imegqwan(5,:))
  if (wproc1) then
    write(151,*)
    write(151,'("Translation limits : ",6I6)')megqwan_tlim(:,1), &
      megqwan_tlim(:,2),megqwan_tlim(:,3)
    call flushifc(151)
  endif
  allocate(idxmegqwan(nwann,nwann,megqwan_tlim(1,1):megqwan_tlim(2,1),&
    megqwan_tlim(1,2):megqwan_tlim(2,2),megqwan_tlim(1,3):megqwan_tlim(2,3)))
  idxmegqwan=-100
  do i=1,nmegqwan
    idxmegqwan(imegqwan(1,i),imegqwan(2,i),imegqwan(3,i),imegqwan(4,i),&
      imegqwan(5,i))=i
  enddo
endif

! setup energy mesh
lr_dw=(lr_w1-lr_w0)/(lr_nw-1)
if (allocated(lr_w)) deallocate(lr_w)
allocate(lr_w(lr_nw))
do i=1,lr_nw
  lr_w(i)=dcmplx(lr_w0+lr_dw*(i-1),lr_eta)/ha2ev
enddo

if (task.eq.401.or.task.eq.402) then
  maxtr_uscrn=1
  ntr_uscrn=(2*maxtr_uscrn+1)**3
  if (allocated(vtl_uscrn)) deallocate(vtl_uscrn)
  allocate(vtl_uscrn(3,ntr_uscrn))
  if (allocated(ivtit_uscrn)) deallocate(ivtit_uscrn)
  allocate(ivtit_uscrn(-maxtr_uscrn:maxtr_uscrn,-maxtr_uscrn:maxtr_uscrn,&
    -maxtr_uscrn:maxtr_uscrn))
  n=0
  do i1=-maxtr_uscrn,maxtr_uscrn
    do i2=-maxtr_uscrn,maxtr_uscrn
      do i3=-maxtr_uscrn,maxtr_uscrn
        n=n+1
        vtl_uscrn(:,n)=(/i1,i2,i3/)
        ivtit_uscrn(i1,i2,i3)=n
      enddo
    enddo
  enddo
  nwloc=mpi_grid_map(lr_nw,dim_k)
  if (allocated(uscrnwan)) deallocate(uscrnwan)
  allocate(uscrnwan(nwann,nwann,ntr_uscrn,nwloc))
  uscrnwan=zzero
  if (allocated(ubarewan)) deallocate(ubarewan)
  allocate(ubarewan(nwann,nwann,ntr_uscrn))
  ubarewan=zzero
endif

! distribute q-vectors along 3-rd dimention
nvq0loc=mpi_grid_map(nvq0,dim_q)

#ifdef _PAPI_
call PAPIF_flops(real_time,cpu_time,fp_ins,mflops,ierr)
#endif

! main loop over q-points
do iqloc=1,nvq0loc
  iq=mpi_grid_map(nvq0,dim_q,loc=iqloc)
  call genmegq(iq)
  if (task.eq.400.or.task.eq.401) call genchi(iq)
  if (task.eq.402) call genubare(iq)
enddo
if (task.eq.401) call write_uscrn
if (task.eq.402) call write_ubare

!
!!-----------------------------------------!
!!    task 400: compute matrix elements    !
!!-----------------------------------------!
!if (task.eq.400) then
!! calculate matrix elements
!  call timer_start(10,reset=.true.)
!  if (wproc1) call timestamp(151,txt='start genmegq')
!  do iq=ivq1,ivq2
!    call genmegq(ivq0m_list(1,iq))
!  enddo
!  call timer_stop(10)
!  if (wproc1) call timestamp(151,txt='stop genmegq')
!  if (wproc1) then
!    write(151,*)
!    write(151,'("Total time for matrix elements : ",F8.2," seconds")')timer_get_value(10)
!    call flushifc(151)
!  endif
!endif
!
!!------------------------------!
!!    task 401: compute chi0    !
!!------------------------------!
!if (task.eq.401) then
!! calculate chi0
!  call timer_start(11,reset=.true.)
!  if (wproc1) call timestamp(151,txt='start genchi0')
!  do iq=ivq1,ivq2
!    call genchi0(ivq0m_list(1,iq))
!  enddo
!  call timer_stop(11)
!  if (wproc1) call timestamp(151,txt='stop genchi0')
!  if (wproc1) then
!    write(151,*)
!    write(151,'("Total time for chi0 : ",F8.2," seconds")')timer_get_value(11)
!    call flushifc(151)
!  endif
!endif
!!---------------------------------------!
!!    task 402: compute me and bare U    !
!!---------------------------------------!
!if (task.eq.402) then
!  if (wproc1) call timestamp(151,txt='start task 402')
!  do iq=ivq1,ivq2
!    call genmegq(ivq0m_list(1,iq))
!    call genubare(ivq0m_list(1,iq))
!  enddo 
!  call write_ubare
!  if (wproc1) call timestamp(151,txt='stop task 402')
!endif
!
!!------------------------------------------!
!!    task 403: compute me, chi0 and chi    !
!!------------------------------------------!
!if (task.eq.403) then
!  if (wproc1) call timestamp(151,txt='start task 403')
!  do iq=ivq1,ivq2
!    call genmegq(ivq0m_list(1,iq))
!    call genchi0(ivq0m_list(1,iq))
!  enddo 
!  if (crpa) call write_uscrn
!  if (wproc1) call timestamp(151,txt='stop task 403')
!endif

#ifdef _PAPI_
call PAPIF_flops(real_time,cpu_time,fp_ins,mflops,ierr)
t1=dble(mflops)
call mpi_grid_reduce(t1)
#endif

if (wproc1) then
  write(151,*)
#ifdef _PAPI_
  write(151,'("Average performance (Gflops/proc) : ",F15.3)')t1/mpi_grid_nproc/1000.d0
#endif  
  write(151,'("Done.")')
  close(151)
endif

deallocate(wfsvmtloc)
deallocate(wfsvitloc)
deallocate(evecfvloc)
deallocate(evecsvloc)
deallocate(ngknr)
deallocate(igkignr)
deallocate(occsvnr)
deallocate(evalsvnr)   
if (wannier_megq) deallocate(wann_c)
return
end
#endif
