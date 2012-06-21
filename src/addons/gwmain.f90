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
integer nvqloc,iqloc,it,ikloc,ik,ist
character*100 qnm,fname
integer nwloc,iwloc,iw
character*8 c1,c2
real(8), allocatable :: vxcnk(:,:)
complex(8) :: dSe, Zf, Eqp

call init0
call init1
if (.not.mpi_grid_in()) return
if (mpi_grid_root()) call timestamp(6,"[gwmain] done init")
call init_q_mesh(8)
call genvq
call genvgq
! create q-directories
!if (mpi_grid_root()) then
!  call system("mkdir -p q")
!  do iq=1,nvq
!    call getqdir(iq,vqm(:,iq),qnm)
!    call system("mkdir -p "//trim(qnm))
!  enddo
!endif
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
  write(151,*)
  do i=1,nvq0
    write(151,'(" vqc : ",3G18.10)')vqc(:,i)
  enddo
  call flushifc(151)
endif
wproc=mpi_grid_root()
! generate wave-functions for entire BZ
call genwfnr(151,tq0bz)
! generate matrix elements <nk|Vxc|nk>
allocate(vxcnk(nstsv,nkptnr))
call genvxcnk(vxcnk)
if (mpi_grid_root()) then
  do ist=1,nstsv
    write(151,'(" ist : ",I4,"    vxcnk : ",G18.10)')ist,vxcnk(ist,1)
  enddo
endif

! setup energy mesh
call gen_w_mesh

! G0W0 calculation
if (gw_mode.eq.0) then
 if (allocated(lr_w)) deallocate(lr_w)
 lr_nw=2
 allocate(lr_w(lr_nw))
 lr_w=0.d0
endif

! distribute frequency points over 1-st dimension
nwloc=mpi_grid_map(lr_nw,dim_k)
! distribute q-vectors along 2-nd dimention
nvqloc=mpi_grid_map(nvq,dim_q)

allocate(gw_self_energy(lr_nw,nstsv,nkptnrloc))
allocate(self_energy_x(nstsv,nkptnrloc))

gw_self_energy=zzero
self_energy_x=zzero

megq_include_bands=chi0_include_bands
! main loop over q-points
do iqloc=1,nvqloc
  iq=mpi_grid_map(nvq,dim_q,loc=iqloc)
  call genmegq(iq,.false.,.true.,.true.)
  call get_adjoint_megqblh(iq)
  call update_self_energy(iq)
  if (mpi_grid_root()) then
    write(151,'("iq : ",I4," out of ",I4)')iq,nvqloc
    call flushifc(151)
  endif
enddo
call mpi_grid_reduce(gw_self_energy(1,1,1),lr_nw*nstsv*nkptnrloc,dims=(/dim_q/))
call mpi_grid_reduce(self_energy_x(1,1),nstsv*nkptnrloc,dims=(/dim_q/))

gw_self_energy(:,:,:)=gw_self_energy(:,:,:)/nkptnr/omega
self_energy_x(:,:)=self_energy_x(:,:)/nkptnr/omega

do iw=1,lr_nw
  gw_self_energy(iw,:,:)=gw_self_energy(iw,:,:)+self_energy_x(:,:) !the total self_energy
enddo

! G0W0 part, modified by I.H. Chu, H.P. Cheng, Jun,2012

if (gw_mode.eq.0) then
  if (mpi_grid_root((/dim_q/))) then
    do ikloc=1,nkptnrloc
      ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
      write(fname,'("Eqp_k",I4.4)')ik
      open(154,file=trim(adjustl(fname)),form="FORMATTED",status="REPLACE")
      do i=1,nstsv
       dSe=gw_self_energy(2,i,ikloc)-gw_self_energy(1,i,ikloc)
       dSe=dSe/dcmplx(del_e)
       Zf=1.d0/(1.d0-dSe)
       Eqp=dcmplx(evalsvnr(i,ik))+ &
           Zf*(gw_self_energy(1,i,ikloc) -dcmplx(vxcnk(i,ik)))
       write(154,'(I5,1X,6f16.6)') i,real(Eqp)*ha2ev, &
                                    evalsvnr(i,ik)*ha2ev,vxcnk(i,ik)*ha2ev, &
                                 real(gw_self_energy(1,i,ikloc))*ha2ev, &
                                 real(self_energy_x(i,ikloc))*ha2ev,real(Zf)

      enddo !i
      write(154,*) "Done!"
      close(154)
    enddo !ikloc
  endif
!!!!!!!!!!!!!!!!!!!!!
else
  if (mpi_grid_root((/dim_q/))) then
    do ikloc=1,nkptnrloc
      ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
      do i=1,nstsv
        write(fname,'("self_energy_k",I4.4,"_b",I4.4,"__.dat")')ik,i
        open(153,file=trim(adjustl(fname)),form="FORMATTED",status="REPLACE")
        do iw=1,lr_nw
          write(153,'(3G18.10)')dreal(lr_w(iw)),dimag(gw_self_energy(iw,i,ikloc)),dreal(gw_self_energy(iw,i,ikloc))
        enddo
        close(153)
      enddo 
    enddo
  endif
endif

if (mpi_grid_root()) then
  call timestamp(151)
  write(151,*) 
  write(151,'("Done.")')
  close(151)
endif
deallocate(lr_w)
deallocate(gw_self_energy)
deallocate(self_energy_x)
deallocate(vxcnk)
return
end
