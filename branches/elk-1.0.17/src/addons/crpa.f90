subroutine crpa
use modmain
use mod_addons_q
use mod_nrkp
implicit none
integer iq,ig,i,n,n1
integer nvqloc,iqloc
real(8) v2(3),vtc(3)
logical lgamma,wproc1
character*100 qnm
complex(8) zt1
integer nwloc

call init0
call init1
if (.not.mpi_grid_in()) return

if (.not.wannier) then
  write(*,*)
  write(*,'("Error(crpa): Wannier functions are off")')
  write(*,*)
  call pstop
endif
wannier_megq=.true.
lgamma=.true.
call init_qbz(lgamma,8)
call init_q_gq
! create q-directories
if (mpi_grid_root()) then
  call system("mkdir -p q")
  do iq=1,nvq
    call getqdir(iq,vqm(:,iq),qnm)
    call system("mkdir -p "//trim(qnm))
  enddo
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
wproc1=.false.
if (mpi_grid_root()) then
  wproc1=.true.
  open(151,file="CRPA.OUT",form="FORMATTED",status="REPLACE")
  call timestamp(151)
  write(151,'("Total number of q-vectors        : ",I6)')nvq
  write(151,'("Total number of processors       : ",I6)')nproc
  write(151,'("MPI grid size                    : ",3I6)')mpi_grid_size
  call flushifc(151)
endif
! generate wave-functions for entire BZ
call genwfnr(151,.true.)
all_wan_ibt=.true.
call getimegqwan(all_wan_ibt)
! setup energy mesh
lr_dw=(lr_w1-lr_w0)/(lr_nw-1)
if (allocated(lr_w)) deallocate(lr_w)
allocate(lr_w(lr_nw))
do i=1,lr_nw
  lr_w(i)=dcmplx(lr_w0+lr_dw*(i-1),lr_eta)/ha2ev
enddo
! distribute frequency points over 1-st dimension
nwloc=mpi_grid_map(lr_nw,dim_k)
! distribute q-vectors along 3-rd dimention
nvqloc=mpi_grid_map(nvq,dim_q)
allocate(uscrnwan(nmegqwan,nwloc))
uscrnwan=zzero

! main loop over q-points
do iqloc=1,nvqloc
  iq=mpi_grid_map(nvq,dim_q,loc=iqloc)
  call genmegq(iq,.true.,.false.)
  call genchi0(iq)
  call genuscrn(iq)
enddo
uscrnwan=uscrnwan/nkptnr/omega
if (wproc1) then
  write(151,*)
  write(151,'("Number of Wannier transitions : ",I6)')nmegqwan
  write(151,'("   n       n''        T                   U_scrn(w=0)")')
  write(151,'(65("-"))')
  do i=1,nmegqwan
    write(151,'(I4,4X,I4,4X,3I3,4X,2G18.10)')imegqwan(:,i),&
      dreal(uscrnwan(i,1)),dimag(uscrnwan(i,1))
  enddo
  close(151)
endif

return
end