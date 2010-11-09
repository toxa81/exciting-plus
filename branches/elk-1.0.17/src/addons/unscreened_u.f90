subroutine unscreened_u
use modmain
use mod_addons_q
use mod_nrkp
implicit none
integer iq,ig,i,n,n1
integer nvqloc,iqloc
real(8) v2(3),vtc(3)
logical wproc1
character*100 qnm
complex(8) zt1

call init0
call init1
if (.not.mpi_grid_in()) return

if (.not.wannier) then
  write(*,*)
  write(*,'("Error(unscreened_u): Wannier functions are off")')
  write(*,*)
  call pstop
endif
wannier_megq=.true.
call init_qbz(tq0bz,1)
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
wproc1=.false.
if (mpi_grid_root()) then
  wproc1=.true.
  open(151,file="UNSCREENED_U.OUT",form="FORMATTED",status="REPLACE")
  call timestamp(151)
  write(151,'("Total number of q-vectors        : ",I6)')nvq
  write(151,'("Total number of processors       : ",I6)')nproc
  write(151,'("MPI grid size                    : ",3I6)')mpi_grid_size
  call flushifc(151)
endif
wproc=wproc1
! generate wave-functions for entire BZ
call genwfnr(151,.false.)
all_wan_ibt=.true.
call getimegqwan(all_wan_ibt)

! distribute q-vectors along 2-nd dimention
nvqloc=mpi_grid_map(nvq,dim_q)
allocate(ubarewan(nmegqwan))
ubarewan=zzero
! main loop over q-points
do iqloc=1,nvqloc
  iq=mpi_grid_map(nvq,dim_q,loc=iqloc)
  call genmegq(iq,.true.,.false.)
  do i=1,nmegqwan
    n=imegqwan(1,i)
    n1=imegqwan(2,i)
    v2=dble(imegqwan(3:5,i))
    call r3mv(avec,v2,vtc)
    zt1=zzero
    do ig=1,ngq(iq)
      zt1=zt1+vhgq(ig,iq)*dconjg(megqwan(idxmegqwan(n,n,0,0,0),ig))*&
        megqwan(idxmegqwan(n1,n1,0,0,0),ig)
    enddo
    ubarewan(i)=ubarewan(i)+zt1*exp(-zi*dot_product(vqc(:,iq),vtc))
  enddo
enddo
ubarewan=ubarewan/nkptnr/omega
if (wproc1) then
  write(151,*)
  write(151,'("Number of Wannier transitions : ",I6)')nmegqwan
  write(151,'("   n       n''        T                     U_bare")')
  write(151,'(65("-"))')
  do i=1,nmegqwan
    write(151,'(I4,4X,I4,4X,3I3,4X,2G18.10)')imegqwan(:,i),&
      dreal(ubarewan(i)),dimag(ubarewan(i))
  enddo
  close(151)
endif

return
end