subroutine gwmain
use modmain
use mod_addons_q
use mod_nrkp
use mod_hdf5
use mod_wannier
use mod_linresp
use mod_expigqr
implicit none
integer iq,i
integer nvqloc,iqloc,it,ikloc,ik
character*100 qnm,fname
integer nwloc,iwloc,iw
character*8 c1,c2

call init0
call init1
if (.not.mpi_grid_in()) return
if (mpi_grid_root()) call timestamp(6,"[gwmain] done init")
call init_qbz(tq0bz,8)
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
! generate radial functions
call genradf
if (mpi_grid_root()) then
  open(151,file="GW.OUT",form="FORMATTED",status="REPLACE")
  call timestamp(151)
  write(151,'("Total number of q-vectors        : ",I6)')nvq
  write(151,'("Total number of processors       : ",I6)')nproc
  write(151,'("MPI grid size                    : ",8I6)')&
    &(mpi_grid_dim_size(i),i=1,mpi_grid_nd)
  call flushifc(151)
endif
wproc=mpi_grid_root()
! generate wave-functions for entire BZ
call genwfnr(151,tq0bz,lmaxvr)
! setup energy mesh
if (lr_nw.eq.1) then
  lr_dw=0.d0
else
  if (timgw) then
    lr_dw=(lr_iw1-lr_iw0)/(lr_nw-1)
  else
    lr_dw=(lr_w1-lr_w0)/(lr_nw-1)
  endif
endif
if (allocated(lr_w)) deallocate(lr_w)
allocate(lr_w(lr_nw))
do i=1,lr_nw
  if (timgw) then
    lr_w(i)=zi*(lr_iw0+lr_dw*(i-1))/ha2ev
  else
    lr_w(i)=dcmplx(lr_w0+lr_dw*(i-1),lr_eta)/ha2ev
  endif
enddo
! distribute frequency points over 1-st dimension
nwloc=mpi_grid_map(lr_nw,dim_k)
! distribute q-vectors along 2-nd dimention
nvqloc=mpi_grid_map(nvq,dim_q)

allocate(self_energy_c(lr_nw,nstsv,nkptnrloc))
self_energy_c=zzero

megq_include_bands=chi0_include_bands
! main loop over q-points
do iqloc=1,nvqloc
  iq=mpi_grid_map(nvq,dim_q,loc=iqloc)
  call genmegq(iq,.true.,.false.)
  call get_adjoint_megqblh
  call update_self_energy(iq)
  if (mpi_grid_root()) then
    write(151,'("iq : ",I4," out of ",I4)')iq,nvqloc
    call flushifc(151)
  endif
enddo
call mpi_grid_reduce(self_energy_c(1,1,1),lr_nw*nstsv*nkptnrloc,dims=(/dim_q/))
self_energy_c(:,:,:)=self_energy_c(:,:,:)/nkptnr/omega

if (mpi_grid_root((/dim_q/))) then
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do i=1,nstsv
      write(fname,'("self_energy_k",I4.4,"_b",I4.4,"__.dat")')ik,i
      open(153,file=trim(adjustl(fname)),form="FORMATTED",status="REPLACE")
      do iw=1,lr_nw
        write(153,'(3G18.10)')dreal(lr_w(iw)),dimag(self_energy_c(iw,i,ikloc)),dreal(self_energy_c(iw,i,ikloc))
      enddo
      close(153)
    enddo 
  enddo
endif

if (mpi_grid_root()) then
  call timestamp(151)
  write(151,*) 
  write(151,'("Done.")')
  close(151)
endif
deallocate(lr_w)
deallocate(self_energy_c)
return
end
