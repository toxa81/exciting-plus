subroutine sic_genvwan
use modmain
use mod_lf
use mod_nrkp
use modxcifc
use mod_addons_q
use mod_hdf5
implicit none

integer n,sz,i,j,n1,ispn,itr,it,itloc
real(8) t1,t2
integer v1l(3)
character*12 c1,c2,c3
character*100 path
real(8) etot_sic

! arrays for Wannier functions
complex(8), allocatable :: wanmt0(:,:,:,:,:)
complex(8), allocatable :: wanir0(:,:,:)

complex(8), allocatable :: vwan_old(:)
logical exist
character*20 fname

complex(8) z1
complex(8), allocatable :: ehart(:),exc(:)

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

!
! Workflow: 
!  compute V_n^{H}(r)   
!  compute W_n(r)
!  compute V_n^{XC}(r)  
!  compute B_n^{XC}(r)
!  compute e_n^{XC}(r)
!  compute E_n = <W_n|V_n^{H}|W_n> + <W_n|e_n^{XC}|W_n>
!  compute V_n = V_n^{H}(r) + V_n^{XC}(r) + \sigma*B_n^{XC}(r)
!  compute W_n*V_n

! initialise universal variables
call init0
call init1
if (.not.mpi_grid_in()) return
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

call lf_init(lf_maxt,dim2)
wproc=mpi_grid_root()
if (wproc) then
  write(fname,'("SIC_",I2.2,".OUT")')iitersic
  open(151,file=trim(adjustl(fname)),form="FORMATTED",status="REPLACE")
endif
if (wproc) then
  sz=lmmaxvr*nrmtmax*natmtot+ngrtot
  sz=16*sz*ntrloc*nspinor*(2*nwann+2)/1024/1024
  write(151,*)
  write(151,'("Required memory for real-space arrays (MB) : ",I6)')sz
  write(151,*)
  call flushifc(151)
endif
! generate wave-functions for all k-points in BZ
call genwfnr(151,.false.)  
! get all Wannier transitions
all_wan_ibt=.true.
call getimegqwan(all_wan_ibt)
! potential of Wannier functions
if (allocated(vwanmt)) deallocate(vwanmt)
allocate(vwanmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor,nwann))
vwanmt=zzero
if (allocated(vwanir)) deallocate(vwanir)
allocate(vwanir(ngrtot,ntrloc,nspinor,nwann))
vwanir=zzero
! generate Hartree potential of Wannier functions
call sic_genvhart
! deallocate unnecessary arrays
deallocate(wfsvmtloc)
deallocate(wfsvitloc)
deallocate(evecfvloc)
deallocate(evecsvloc)
deallocate(wann_c)
! restore wproc
wproc=mpi_grid_root()
if (wproc) then
  write(151,*)
  write(151,'("ngqmax : ",I4)')ngqmax
  write(151,'("time for q-vectors : ",F8.3)')timer_get_value(10)
  write(151,'("time for Hartree potential : ",F8.3)')timer_get_value(11)
endif
! both components of spinor wave-function feel the same Coulomb potential
if (spinpol) then
  do n=1,nwann
    vwanmt(:,:,:,:,2,n)=vwanmt(:,:,:,:,1,n)
    vwanir(:,:,2,n)=vwanir(:,:,1,n)
  enddo
endif
! Wannier functions
if (allocated(wanmt)) deallocate(wanmt)
allocate(wanmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor,nwann))
if (allocated(wanir)) deallocate(wanir)
allocate(wanir(ngrtot,ntrloc,nspinor,nwann))
! generate Wannier functions on a mesh
call timer_reset(1)
call timer_reset(2)
allocate(wanmt0(lmmaxvr,nrmtmax,natmtot,nspinor,nwann))
allocate(wanir0(ngrtot,nspinor,nwann))
do itloc=1,ntrloc
  itr=mpi_grid_map(ntr,dim_t,loc=itloc)
  call sic_genwann(vtl(1,itr),ngknr,vgkcnr,igkignr,wanmt0,wanir0)
  do ispn=1,nspinor
    do n=1,nwann
      wanmt(:,:,:,itloc,ispn,n)=wanmt0(:,:,:,ispn,n)
      wanir(:,itloc,ispn,n)=wanir0(:,ispn,n)
    enddo !n
  enddo !ispn
enddo !itr
deallocate(wanmt0,wanir0)
if (wproc) then
  write(151,*)
  write(151,'("Wann MT part : ",F8.3)')timer_get_value(1)
  write(151,'("Wann IT part : ",F8.3)')timer_get_value(2)
  call flushifc(151)
endif
! check orthonormality
t1=0.d0
t2=0.d0
do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  v1l(:)=imegqwan(3:5,i)
  z1=0
  do ispn=1,nspinor
    z1=z1+lf_dotlf(.true.,v1l,wanmt(1,1,1,1,ispn,n),wanir(1,1,ispn,n),&
      wanmt(1,1,1,1,ispn,n1),wanir(1,1,ispn,n1))
  enddo
  if (n.eq.n1.and.v1l(1).eq.0.and.v1l(2).eq.0.and.v1l(3).eq.0) then
    z1=z1-zone
  endif
  t2=max(t2,abs(z1))
  t1=t1+abs(z1)
enddo
if (wproc) then
  write(151,*)
  write(151,'("Maximum deviation from norm : ",F12.6)')t2
  write(151,'("Average deviation from norm : ",F12.6)')t1/nmegqwan
  call flushifc(151)
endif
if (wproc) then
  call timestamp(151,'done with Wannier functions')
endif
! compute product V^{H}(r)*W(r)
do n=1,nwann
  do ispn=1,nspinor
    call lf_prod(zone,vwanmt(1,1,1,1,ispn,n),vwanir(1,1,ispn,n),&
      wanmt(1,1,1,1,ispn,n),wanir(1,1,ispn,n),zzero,&
      vwanmt(1,1,1,1,ispn,n),vwanir(1,1,ispn,n))
  enddo
enddo
! compute <W_n|V^H|W_n>
allocate(ehart(nwann))
ehart=zzero
do n=1,nwann
  do ispn=1,nspinor
    ehart(n)=ehart(n)+lf_dotlf(.true.,(/0,0,0/),vwanmt(1,1,1,1,ispn,n),&
      vwanir(1,1,ispn,n),wanmt(1,1,1,1,ispn,n),wanir(1,1,ispn,n))
  enddo
enddo
! XC part
allocate(exc(nwann))
exc=zzero
call timer_start(12,reset=.true.)
call sic_genvxc(exc)
call timer_stop(12)
if (wproc) then
  write(151,'("time for XC potential : ",F8.3)')timer_get_value(12)
endif

etot_sic=0.d0
do n=1,nwann
  etot_sic=etot_sic+dreal(ehart(n))+dreal(exc(n))
enddo
if (wproc) then
  write(151,'(2X,"wann",16X,"E_n^{H}",33X,"E_n^{XC}")')
  write(151,'(84("-"))')
  do n=1,nwann
    write(151,'(I4,4X,2G18.10,4X,2G18.10)')n,dreal(ehart(n)),dimag(ehart(n)),&
      dreal(exc(n)),dimag(exc(n))
  enddo
  write(151,*)
  write(151,'("Total energy correction : ",G18.10)')etot_sic
  call flushifc(151)
endif

if (allocated(vwan)) deallocate(vwan)
allocate(vwan(nmegqwan))
vwan=zzero
! compute matrix elements of SIC potential
! vwan = <w_n|v_n|w_{n1,T}>
do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  v1l(:)=imegqwan(3:5,i)
  do ispn=1,nspinor
    vwan(i)=vwan(i)+lf_dotlf(.true.,v1l,vwanmt(1,1,1,1,ispn,n),&
      vwanir(1,1,ispn,n),wanmt(1,1,1,1,ispn,n1),wanir(1,1,ispn,n1))
  enddo
enddo
if (wproc) then
  call timestamp(151,'done with matrix elements')
endif
if (wproc) then
  write(151,*)
  write(151,'("Number of Wannier transitions : ",I6)')nmegqwan
  write(151,'("Matrix elements of SIC potential &
    &(n n1  T  <w_n|v_n|w_{n1,T}>)")')
  do i=1,nmegqwan
    write(151,'(I4,4X,I4,4X,3I3,4X,2G18.10)')imegqwan(:,i),&
      dreal(vwan(i)),dimag(vwan(i))
  enddo
endif

t2=0.d0
do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  v1l(:)=imegqwan(3:5,i)
  j=idxmegqwan(n1,n,-v1l(1),-v1l(2),-v1l(3))
  t2=max(t2,abs(vwan(i)-dconjg(vwan(j))))
enddo
if (wproc) then
  write(151,*)
  write(151,'("Maximum deviation from ""localization criterion"" : ",F12.6)')t2
  write(151,*)
  write(151,'("Diagonal matrix elements")')
  write(151,'(2X,"wann",18X,"V_n")')
  write(151,'(44("-"))')
  do n=1,nwann
    j=idxmegqwan(n,n,0,0,0)
    write(151,'(I4,4X,2G18.10)')n,dreal(vwan(j)),dimag(vwan(j))
  enddo  
  call flushifc(151)
endif
if (wproc) then
  inquire(file="sic.hdf5",exist=exist)
  if (exist) then
    allocate(vwan_old(nmegqwan))
    call hdf5_read("sic.hdf5","/","vwan",vwan_old(1),(/nmegqwan/))
    t1=0.d0
    do i=1,nmegqwan
      t1=t1+abs(vwan(i)-vwan_old(i))**2
    enddo
    t1=sqrt(t1/nmegqwan)
    write(151,*)
    write(151,'("SIC matrix elements RMS difference :",G18.10)')t1
  endif
endif

if (wproc) then
  call hdf5_create_file("sic.hdf5")
  call hdf5_create_group("sic.hdf5","/","wann")
  do n=1,nwann
    path="/wann"
    write(c1,'("n",I4.4)')n
    call hdf5_create_group("sic.hdf5",path,trim(adjustl(c1)))   
    path=trim(path)//"/"//trim(adjustl(c1))
    do ispn=1,nspinor
      write(c3,'("s",I4.4)')ispn
      call hdf5_create_group("sic.hdf5",path,trim(adjustl(c3)))   
      path=trim(path)//"/"//trim(adjustl(c3))      
      do it=1,ntr
        write(c2,'("t",I4.4)')it
        call hdf5_create_group("sic.hdf5",path,trim(adjustl(c2)))
      enddo
    enddo
  enddo
  call hdf5_write("sic.hdf5","/","nmegqwan",nmegqwan)
  call hdf5_write("sic.hdf5","/","imegqwan",imegqwan(1,1),(/5,nmegqwan/))
  call hdf5_write("sic.hdf5","/","vwan",vwan(1),(/nmegqwan/))
  !call hdf5_write("sic.hdf5","/","hwan",hwan(1),(/nmegqwan/))
endif
if (mpi_grid_side(dims=(/dim_t/))) then
  do i=0,mpi_grid_size(dim_t)-1
    if (mpi_grid_x(dim_t).eq.i) then
      do itloc=1,ntrloc
        itr=mpi_grid_map(ntr,dim_t,loc=itloc)
        do ispn=1,nspinor
          do n=1,nwann
            write(c1,'("n",I4.4)')n
            write(c2,'("t",I4.4)')itr
            write(c3,'("s",I4.4)')ispn
            path="/wann/"//trim(adjustl(c1))//"/"//trim(adjustl(c3))//"/"//&
              trim(adjustl(c2))
            call hdf5_write("sic.hdf5",path,"vwanmt",&
              vwanmt(1,1,1,itloc,ispn,n),(/lmmaxvr,nrmtmax,natmtot/))
            call hdf5_write("sic.hdf5",path,"vwanir",&
              vwanir(1,itloc,ispn,n),(/ngrtot/))
!            call hdf5_write("sic.hdf5",path,"wanmt",&
!              wanmt(1,1,1,itloc,ispn,n),(/lmmaxvr,nrmtmax,natmtot/))
!            call hdf5_write("sic.hdf5",path,"wanir",&
!              wanir(1,itloc,ispn,n),(/ngrtot/))
          enddo
        enddo
      enddo
    endif
    call mpi_grid_barrier(dims=(/dim_t/))
  enddo
endif
if (wproc) close(151)
call bstop
deallocate(vwan)
deallocate(wanmt,wanir,vwanmt,vwanir)
return
end

