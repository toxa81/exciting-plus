#ifdef _HDF5_
subroutine response
use modmain
use hdf5
use mod_nrkp
use mod_addons_q
implicit none
#ifdef _PAPI_
integer ierr
include 'f90papi.h'
real real_time,cpu_time,mflops
integer*8 fp_ins
#endif
integer i,j,iq
integer nvqloc,iqloc
character*100 qnm
logical wproc1,lpmat
logical, external :: wann_diel
! MPI grid:
!   (matrix elements) : (1) k-points x (2) q-points x 
!                     x (3) interband transitions 
!   (response) : (1) energy mesh x (2) q-points x (3) number of fxc kernels
!
if (lrtype.eq.1.and..not.spinpol) then
  write(*,*)
  write(*,'("Error(response): can''t do magnetic response for unpolarized &
    &ground state")')
  write(*,*)
  call pstop
endif
if (nvq.eq.0) then
  write(*,*)
  write(*,'("Error(response): no q-vectors given")')
  write(*,*)
  call pstop
endif   
if (.not.wannier) wannier_chi0_chi=.false.
if (.not.spinpol) megqwan_afm=.false.
if (wannier_chi0_chi) wannier_megq=.true.

! this is enough for matrix elements
!lmaxvr=5

! initialise universal variables
call init0
call init1
if (.not.mpi_grid_in()) then
  write(*,'("Info(response) processor ",I6," is not in grid")')iproc
  return
endif

! check if momentum matrix is required
lpmat=.false.
do j=1,nvq
  if (all(vqm(:,j).eq.0)) lpmat=.true.
enddo
! create q-directories
if (mpi_grid_root()) then
  call system("mkdir -p q")
  do iq=1,nvq
    call getqdir(iq,vqm(:,iq),qnm)
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
  write(151,'("Total number of q-vectors        : ",I6)')nvq
  write(151,'("Total number of processors       : ",I6)')nproc
  write(151,'("MPI grid size                    : ",3I6)')mpi_grid_size
  write(151,'("Wannier functions                : ",L1)')wannier
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
call getufr
call genufrp
! read Fermi energy
if (mpi_grid_root()) call readfermi
call mpi_grid_bcast(efermi)
! generate wave-functions for entire BZ
call genwfnr(151,lpmat)

if (wannier_megq) then
  all_wan_ibt=.false.
  call getimegqwan(all_wan_ibt)
  if (wproc1) then
    write(151,*)
    write(151,'("Number of Wannier transitions : ",I6)')nmegqwan
    write(151,'("List of Wannier transitions (n n1 T) ")')
    do i=1,nmegqwan
      write(151,'(I4,4X,I4,4X,3I3)')imegqwan(:,i)
    enddo
    call timestamp(151)
    write(151,*)
    write(151,'("Translation limits : ",6I6)')megqwan_tlim(:,1), &
      megqwan_tlim(:,2),megqwan_tlim(:,3)
    call flushifc(151)
  endif
endif
! setup energy mesh
lr_dw=(lr_w1-lr_w0)/(lr_nw-1)
if (allocated(lr_w)) deallocate(lr_w)
allocate(lr_w(lr_nw))
do i=1,lr_nw
  lr_w(i)=dcmplx(lr_w0+lr_dw*(i-1),lr_eta)/ha2ev
enddo
call init_q_gq
! distribute q-vectors along 3-rd dimention
nvqloc=mpi_grid_map(nvq,dim_q)
#ifdef _PAPI_
call PAPIF_flops(real_time,cpu_time,fp_ins,mflops,ierr)
#endif
do iqloc=1,nvqloc
  iq=mpi_grid_map(nvq,dim_q,loc=iqloc)
  call genmegq(iq,.true.,.true.)
  call genchi0(iq)
  call genchi(iq)
enddo
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
