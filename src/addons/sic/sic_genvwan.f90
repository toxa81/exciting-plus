subroutine sic_genvwan
use modmain
use mod_lf
use mod_nrkp
use modxcifc
use mod_addons_q
use mod_hdf5
implicit none
integer n,sz,i,j,n1,ispn,itr,itloc
real(8) t1,t2
integer v1l(3)
!complex(8), allocatable :: vhwanmt(:,:,:,:,:)
!complex(8), allocatable :: vhwanir(:,:,:)
! potential (Hartree+XC) of Wannier function charge density
real(8), allocatable :: vhxcmt(:,:,:,:,:)
real(8), allocatable :: vhxcir(:,:,:)
! Wannier functions
complex(8), allocatable :: wanmt(:,:,:,:,:,:)
complex(8), allocatable :: wanir(:,:,:,:)
complex(8), allocatable :: wanmt0(:,:,:,:,:)
complex(8), allocatable :: wanir0(:,:,:)
complex(8), allocatable :: vwan_old(:)
complex(8), allocatable :: ehart(:),exc(:)
complex(8) z1
logical exist

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

wproc=mpi_grid_root()
if (wproc) then
  open(151,file="SIC.OUT",form="FORMATTED",status="REPLACE")
endif
if (wproc) then
  sz=lmmaxvr*nrmtmax*natmtot+ngrtot
  sz=16*sz*ntrloc*nspinor*(2*nwann+2)/1024/1024
  write(151,*)
  write(151,'("Required memory for real-space arrays (MB) : ",I6)')sz
  write(151,*)
  write(151,'("cutoff radius for WF : ",F12.6)')wann_r_cutoff
  write(151,'("number of translations : ",I4)')ntr
  do i=1,ntr
    write(151,'("  i : ",I4,"    vtl(i) : ",3I4)')i,vtl(:,i)
  enddo
  call flushifc(151)
endif

! generate wave-functions for all k-points in BZ
call genwfnr(151,.false.)  
! get all Wannier transitions
all_wan_ibt=.true.
call getimegqwan(all_wan_ibt)
wvmt=zzero
wvir=zzero
!-------------------!
! Hartree potential !
!-------------------!
! generate Hartree potential of Wannier functions
! use wvmt,wvir arrays as temporary
call sic_genvhart(wvmt,wvir)
! deallocate unnecessary arrays
deallocate(wfsvmtnrloc)
deallocate(wfsvitnrloc)
deallocate(wanncnrloc)
! restore wproc
wproc=mpi_grid_root()
if (wproc) then
  write(151,*)
  write(151,'("ngqmax : ",I4)')ngqmax
  write(151,'("time for q-vectors : ",F8.3)')timer_get_value(10)
  write(151,'("time for Hartree potential : ",F8.3)')timer_get_value(11)
  write(151,'("average imaginary part (mt,ir) : ",2G18.10)') &
    sum(abs(dimag(wvmt)))/lmmaxvr/nrmtmax/natmtot/ntrloc/nwann,&
    sum(abs(dimag(wvir)))/ngrtot/ntrloc/nwann
endif
allocate(vhxcmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nwann))
allocate(vhxcir(ngrtot,ntrloc,nwann))
do n=1,nwann
  vhxcmt(:,:,:,:,n)=dreal(wvmt(:,:,:,:,1,n))
  vhxcir(:,:,n)=dreal(wvir(:,:,1,n))
enddo
!-------------------!
! Wannier functions !
!-------------------!
allocate(wanmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor,nwann))
allocate(wanir(ngrtot,ntrloc,nspinor,nwann))
! generate Wannier functions on a mesh
call timer_reset(1)
call timer_reset(2)
allocate(wanmt0(lmmaxvr,nrmtmax,natmtot,nspinor,nwann))
allocate(wanir0(ngrtot,nspinor,nwann))
do itloc=1,ntrloc
  itr=mpi_grid_map(ntr,dim_t,loc=itloc)
  call sic_genwann(vtl(1,itr),ngknr,igkignr,wanmt0,wanir0)
  do ispn=1,nspinor
    do n=1,nwann
      wanmt(:,:,:,itloc,ispn,n)=wanmt0(:,:,:,ispn,n)
      wanir(:,itloc,ispn,n)=wanir0(:,ispn,n)
    enddo !n
  enddo !ispn
enddo !itr
deallocate(wanmt0,wanir0)
deallocate(wann_unkmt)
deallocate(wann_unkit)
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
    z1=z1+lf_dot_ll(.true.,v1l,wanmt(1,1,1,1,ispn,n),wanir(1,1,ispn,n),&
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
!------------------------------!
! Hartree energy <W_n|V^H|W_n> !
!------------------------------!
allocate(ehart(nwann))
ehart=zzero
do n=1,nwann
  do ispn=1,nspinor
    ehart(n)=ehart(n)+lf_intgr_zdz(wanmt(1,1,1,1,ispn,n),wanir(1,1,ispn,n),&
      vhxcmt(1,1,1,1,n),vhxcir(1,1,n),(/0,0,0/),wanmt(1,1,1,1,ispn,n),&
      wanir(1,1,ispn,n))
  enddo
  write(*,*)ehart(n)
enddo
call pstop
!!-------------------------!
!! XC potential and energy !
!!-------------------------!
!allocate(exc(nwann))
!exc=zzero
!call timer_start(12,reset=.true.)
!call sic_genvxc(exc)
!call timer_stop(12)
!if (wproc) then
!  write(151,'("time for XC potential : ",F8.3)')timer_get_value(12)
!endif
!sic_etot_correction=0.d0
!ehart(:)=-1.d0*ehart(:)
!exc(:)=-1.d0*exc(:)
!do n=1,nwann
!  sic_etot_correction=sic_etot_correction+dreal(ehart(n))+dreal(exc(n))
!enddo
!if (wproc) then
!  write(151,'(2X,"wann",16X,"E_n^{H}",33X,"E_n^{XC}")')
!  write(151,'(84("-"))')
!  do n=1,nwann
!    write(151,'(I4,4X,2G18.10,4X,2G18.10)')n,dreal(ehart(n)),dimag(ehart(n)),&
!      dreal(exc(n)),dimag(exc(n))
!  enddo
!  write(151,*)
!  write(151,'("Total energy correction : ",G18.10)')sic_etot_correction
!  call flushifc(151)
!endif
!!----------------------------------!
!! matrix elements of SIC potential !
!!----------------------------------!
!do n=1,nwann
!  do ispn=1,nspinor
!    call lf_mult_zd(-zone,wanmt(1,1,1,1,ispn,n),wanir(1,1,ispn,n), &
!      vwanmt(1,1,1,1,ispn,n),vwanir(1,1,ispn,n),vwanmt_(1,1,1,1,ispn,n),&
!      vwanir_(1,1,ispn,n))    
!  enddo
!enddo
!if (allocated(vwan)) deallocate(vwan)
!allocate(vwan(nmegqwan))
!vwan=zzero
!! compute matrix elements of SIC potential
!! vwan = <w_n|v_n|w_{n1,T}>
!do i=1,nmegqwan
!  n=imegqwan(1,i)
!  n1=imegqwan(2,i)
!  v1l(:)=imegqwan(3:5,i)
!  do ispn=1,nspinor    
!    vwan(i)=vwan(i)+lf_intgr_zz(vwanmt_(1,1,1,1,ispn,n),vwanir_(1,1,ispn,n),&
!      v1l,wanmt(1,1,1,1,ispn,n1),wanir(1,1,ispn,n1))
!  enddo
!enddo
!if (wproc) then
!  call timestamp(151,'done with matrix elements')
!endif
!if (wproc) then
!  write(151,*)
!  write(151,'("Number of Wannier transitions : ",I6)')nmegqwan
!  write(151,'("Matrix elements of SIC potential &
!    &(n n1  T  <w_n|v_n|w_{n1,T}>)")')
!  do i=1,nmegqwan
!    write(151,'(I4,4X,I4,4X,3I3,4X,2G18.10)')imegqwan(:,i),&
!      dreal(vwan(i)),dimag(vwan(i))
!  enddo
!endif
!
!t1=0.d0
!t2=0.d0
!do i=1,nmegqwan
!  n=imegqwan(1,i)
!  n1=imegqwan(2,i)
!  v1l(:)=imegqwan(3:5,i)
!  j=idxmegqwan(n1,n,-v1l(1),-v1l(2),-v1l(3))
!  t1=t1+abs(vwan(i)-dconjg(vwan(j)))
!  t2=max(t2,abs(vwan(i)-dconjg(vwan(j))))
!enddo
!if (wproc) then
!  write(151,*)
!  write(151,'("Maximum deviation from ""localization criterion"" : ",F12.6)')t2
!  write(151,'("Average deviation from ""localization criterion"" : ",F12.6)')t1/nmegqwan
!  write(151,*)
!  write(151,'("Diagonal matrix elements")')
!  write(151,'(2X,"wann",18X,"V_n")')
!  write(151,'(44("-"))')
!  do n=1,nwann
!    j=idxmegqwan(n,n,0,0,0)
!    write(151,'(I4,4X,2G18.10)')n,dreal(vwan(j)),dimag(vwan(j))
!  enddo  
!  call flushifc(151)
!endif
!if (wproc) then
!  inquire(file="sic.hdf5",exist=exist)
!  if (exist) then
!    allocate(vwan_old(nmegqwan))
!    call hdf5_read("sic.hdf5","/","vwan",vwan_old(1),(/nmegqwan/))
!    t1=0.d0
!    do i=1,nmegqwan
!      t1=t1+abs(vwan(i)-vwan_old(i))**2
!    enddo
!    t1=sqrt(t1/nmegqwan)
!    write(151,*)
!    write(151,'("SIC matrix elements RMS difference :",G18.10)')t1
!    deallocate(vwan_old)
!  endif
!endif
!call sic_writevwan
!if (wproc) close(151)
!deallocate(vwan)
!deallocate(ehart)
!deallocate(exc)
!deallocate(vwanmt,vwanir)
!deallocate(wanmt,wanir)
return
end

