subroutine crpa
use modmain
use mod_addons_q
use mod_nrkp
use mod_hdf5
use mod_wannier
use mod_linresp
use mod_expigqr
implicit none
integer iq,i
integer nvqloc,iqloc,ntrloc,it,itloc
logical wproc1
character*100 qnm,fu4
integer nwloc,iwloc,iw
character*8 c1,c2
logical exist
integer*8, allocatable :: hw_values(:)
logical, parameter :: u4scrn=.true.

call init0
call init1
if (.not.mpi_grid_in()) return
if (mpi_grid_root()) call timestamp(6,"[crpa] done init")
if (.not.wannier) then
  write(*,*)
  write(*,'("Error(crpa): Wannier functions are off")')
  write(*,*)
  call pstop
endif
call papi_timer_start(pt_crpa_tot1)
wannier_megq=.true.
if (u4scrn) then
  call init_qbz(tq0bz,8)
else
  call init_qbz(tq0bz,1)
endif
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
  open(151,file="CRPA.OUT",form="FORMATTED",status="REPLACE")
  call timestamp(151)
  write(151,'("Total number of q-vectors        : ",I6)')nvq
  write(151,'("Total number of processors       : ",I6)')nproc
  write(151,'("MPI grid size                    : ",3I6)')mpi_grid_size
  call flushifc(151)
endif
wproc=wproc1
inquire(file="wfnrkp.hdf5",exist=exist)
if (exist) then
  call timer_start(1,reset=.true.)
  call drc_read_wf
  call timer_stop(1)
  if (wproc1) then
    write(151,*)
    write(151,'("drc_read_wf done in ",F8.2," seconds")')timer_get_value(1)
    call flushifc(151)
  endif
else
! generate wave-functions for entire BZ
  call genwfnr(151,tq0bz)
endif
call genwantran(megqwantran,megqwan_mindist,megqwan_maxdist,allwt=.true.)
! setup energy mesh
if (.not.u4scrn) lr_nw=1
if (lr_nw.eq.1) then
  lr_dw=0.d0
else
  lr_dw=(lr_w1-lr_w0)/(lr_nw-1)
endif
if (allocated(lr_w)) deallocate(lr_w)
allocate(lr_w(lr_nw))
do i=1,lr_nw
  lr_w(i)=dcmplx(lr_w0+lr_dw*(i-1),lr_eta)/ha2ev
enddo
! distribute frequency points over 1-st dimension
nwloc=mpi_grid_map(lr_nw,dim_k)
! distribute q-vectors along 2-nd dimention
nvqloc=mpi_grid_map(nvq,dim_q)
! distribute translations along 3-nd dimention
ntrloc=mpi_grid_map(megqwantran%ntr,dim_b)

if (mpi_grid_root()) then
  write(151,*)
  write(151,'("Number of Wannier transitions : ",I6)')megqwantran%nwt
  write(151,'("Number of Wannier translations : ",I6)')megqwantran%ntr
  write(*,*)
  write(*,'("[crpa] size of 4-index U matrix : ",I6," Mb")') &
        int(16.d0*megqwantran%nwt*megqwantran%nwt*ntrloc*nwloc/1048576.d0)
endif
call mpi_grid_barrier()
allocate(u4(megqwantran%nwt,megqwantran%nwt,ntrloc,nwloc))
u4=zzero
call papi_timer_start(pt_crpa_tot2)
! main loop over q-points
do iqloc=1,nvqloc
  iq=mpi_grid_map(nvq,dim_q,loc=iqloc)
  call genmegq(iq,.true.,.false.)
  call genu4(iq,nwloc, ntrloc,u4scrn)
enddo
do iwloc=1,nwloc
  do itloc=1,ntrloc
    call mpi_grid_reduce(u4(1,1,itloc,iwloc),megqwantran%nwt*megqwantran%nwt,&
      dims=(/dim_q/))
  enddo
enddo
call papi_timer_stop(pt_crpa_tot2)
call papi_timer_stop(pt_crpa_tot1)

allocate(hw_values(0:papi_ncounters))
call papi_timer_read(pt_crpa_tot1,hw_values)
call mpi_grid_reduce(hw_values(0),1+papi_ncounters)
if (wproc1) call papi_report(151,hw_values,"pt_crpa_tot1")
call papi_timer_read(pt_crpa_tot2,hw_values)
call mpi_grid_reduce(hw_values(0),1+papi_ncounters)
if (wproc1) call papi_report(151,hw_values,"pt_crpa_tot2")
call papi_timer_read(pt_megq,hw_values)
call mpi_grid_reduce(hw_values(0),1+papi_ncounters)
if (wproc1) call papi_report(151,hw_values,"pt_megq")
call papi_timer_read(pt_megqblh,hw_values)
call mpi_grid_reduce(hw_values(0),1+papi_ncounters)
if (wproc1) call papi_report(151,hw_values,"pt_megqblh")
call papi_timer_read(pt_megqblh_mt,hw_values)
call mpi_grid_reduce(hw_values(0),1+papi_ncounters)
if (wproc1) call papi_report(151,hw_values,"pt_megqblh_mt")
call papi_timer_read(pt_megqblh_it,hw_values)
call mpi_grid_reduce(hw_values(0),1+papi_ncounters)
if (wproc1) call papi_report(151,hw_values,"pt_megqblh_it")
call papi_timer_read(pt_megqblh2,hw_values)
call mpi_grid_reduce(hw_values(0),1+papi_ncounters)
if (wproc1) call papi_report(151,hw_values,"pt_megqblh2")
call papi_timer_read(pt_chi0_zgemm,hw_values)
call mpi_grid_reduce(hw_values(0),1+papi_ncounters)
if (wproc1) call papi_report(151,hw_values,"pt_chi0_zgemm")
call papi_timer_read(pt_chi0,hw_values)
call mpi_grid_reduce(hw_values(0),1+papi_ncounters)
if (wproc1) call papi_report(151,hw_values,"pt_chi0")
call papi_timer_read(pt_uscrn,hw_values)
call mpi_grid_reduce(hw_values(0),1+papi_ncounters)
if (wproc1) call papi_report(151,hw_values,"pt_uscrn")
call papi_timer_read(pt_vscrn,hw_values)
call mpi_grid_reduce(hw_values(0),1+papi_ncounters)
if (wproc1) call papi_report(151,hw_values,"pt_vscrn")
deallocate(hw_values)

if (mpi_grid_side(dims=(/dim_k,dim_b/)).and.nwloc.gt.0) then
  write(fu4,'("u4_",I4.4,".hdf5")')mpi_grid_x(dim_k)
  if (mpi_grid_root(dims=(/dim_b/))) then
    call hdf5_create_file(trim(fu4))
    call hdf5_create_group(trim(fu4),"/","iwloc")
    call hdf5_create_group(trim(fu4),"/","parameters")
    call hdf5_write(fu4,"/parameters","nwantot",nwantot)  
    call hdf5_write(fu4,"/parameters","nw",lr_nw)
    call hdf5_write(fu4,"/parameters","nwloc",nwloc)
    call hdf5_write(fu4,"/parameters","x",mpi_grid_x(dim_k))
    call hdf5_write(fu4,"/parameters","size",mpi_grid_size(dim_k))
    call hdf5_write(fu4,"/parameters","nwt",megqwantran%nwt)
    call hdf5_write(fu4,"/parameters","iwt",megqwantran%iwt(1,1),&
      (/5,megqwantran%nwt/))  
    call hdf5_write(fu4,"/parameters","ngq",ngq(1))
    call hdf5_write(fu4,"/parameters","ngridk",ngridk(1),(/3/))
    call hdf5_write(fu4,"/parameters","ntr",megqwantran%ntr)
    call hdf5_write(fu4,"/parameters","vtr",megqwantran%vtr(1,1),&
      (/3,megqwantran%ntr/))
  endif
  do iwloc=1,nwloc
    iw=mpi_grid_map(lr_nw,dim_k,loc=iwloc)
    write(c1,'(I8.8)')iwloc
    if (mpi_grid_root(dims=(/dim_b/))) then
      call hdf5_create_group(trim(fu4),"/iwloc",c1)
      do it=1,megqwantran%ntr
        write(c2,'("t",I7.7)')it
        call hdf5_create_group(trim(fu4),"/iwloc/"//c1,c2)
      enddo
      call hdf5_write(fu4,"/iwloc/"//c1,"iw",iw)
      call hdf5_write(fu4,"/iwloc/"//c1,"w",dreal(lr_w(iw)))
    endif
    do i=0,mpi_grid_size(dim_b)-1
      if (i.eq.mpi_grid_x(dim_b)) then
        do itloc=1,ntrloc
          it=mpi_grid_map(megqwantran%ntr,dim_b,loc=itloc)
          write(c2,'("t",I7.7)')it
          call hdf5_write(fu4,"/iwloc/"//c1//"/"//c2,"u4",u4(1,1,itloc,iwloc),&
            (/megqwantran%nwt,megqwantran%nwt/))
        enddo
      endif
      call mpi_grid_barrier(dims=(/dim_b/))
    enddo
  enddo
endif
if (wproc1) then
!  write(151,*)
!  write(151,'("Number of Wannier transitions : ",I6)')megqwantran%nwt
!  write(151,'("Number of Wannier translations : ",I6)')megqwantran%ntr
!  write(151,*)
!  write(151,'("n   n''   Re,Im U_{n,n''T=0} values for w=",2G18.10)')dreal(lr_w(1)),&
!    dimag(lr_w(1))
!  do i=1,megqwantran%nwan
!    n=megqwantran%iwan(i)
!    j=megqwantran%iwtidx(n,n,0,0,0)
!    do i1=1,megqwantran%nwan
!      n1=megqwantran%iwan(i1)
!      j1=megqwantran%iwtidx(n1,n1,0,0,0)
!      if (j.ne.0.and.j1.ne.0) then
!        write(151,'(2I4,2G18.10)')n,n1,dreal(u4(
!      else
!        write(151,'("wrong index")')
!        write(151,'(" j, j1 : ",2I6)')j,j1
!      endif
!        
!  enddo
!  write(151,'("screened U_{n,n''T}(w=0)")')    
!  write(151,'(65("-"))')
!    call printwanntrans(151,uscrnwan(1,1))
!  call timestamp(151)
!  write(151,*)
!  write(151,'("screened J_{n,n''T}(w=0)")')    
!  write(151,'(65("-"))')
!    call printwanntrans(151,jscrnwan(1,1))
  call timestamp(151)
  write(151,*) 
  write(151,'("Done.")')
  close(151)
endif
deallocate(lr_w)
deallocate(u4)
return
end