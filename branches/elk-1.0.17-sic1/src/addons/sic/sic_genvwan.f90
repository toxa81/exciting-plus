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
! potential (Hartree+XC) of Wannier function charge density
real(8), allocatable :: vhxcmt(:,:,:,:,:,:)
real(8), allocatable :: vhxcir(:,:,:,:)
! Wannier functions
complex(8), allocatable :: wanmt(:,:,:,:,:,:)
complex(8), allocatable :: wanir(:,:,:,:)
complex(8), allocatable :: wanmt0(:,:,:,:,:)
complex(8), allocatable :: wanir0(:,:,:)
complex(8), allocatable :: vwanme_old(:)
complex(8), allocatable :: ehart(:),exc(:)
complex(8) z1
complex(8), allocatable :: norm(:)
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
  call timestamp(151,'done with Hartree potential')
endif
allocate(vhxcmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor,nwann))
allocate(vhxcir(ngrtot,ntrloc,nspinor,nwann))
do ispn=1,nspinor
  do n=1,nwann
    vhxcmt(:,:,:,:,ispn,n)=dreal(wvmt(:,:,:,:,1,n))
    vhxcir(:,:,ispn,n)=dreal(wvir(:,:,1,n))
  enddo
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
allocate(norm(nwann))
! check orthonormality
t1=0.d0
t2=0.d0
do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  v1l(:)=imegqwan(3:5,i)
  z1=0
  do ispn=1,nspinor
    z1=z1+lf_dot_lf(.true.,wanmt(1,1,1,1,ispn,n),wanir(1,1,ispn,n),v1l,&
      wanmt(1,1,1,1,ispn,n1),wanir(1,1,ispn,n1))
  enddo
  if (n.eq.n1.and.v1l(1).eq.0.and.v1l(2).eq.0.and.v1l(3).eq.0) then
    norm(n)=z1
    z1=z1-zone
  endif
  t2=max(t2,abs(z1))
  t1=t1+abs(z1)
enddo
if (wproc) then
  write(151,*)
  do n=1,nwann
    write(151,'("n : ",I4,"   norm : ",2G18.10)')n,dreal(norm(n)),&
      dimag(norm(n))
  enddo
  write(151,*)
  write(151,'("Maximum deviation from norm : ",F12.6)')t2
  write(151,'("Average deviation from norm : ",F12.6)')t1/nmegqwan
  call flushifc(151)
endif
if (wproc) then
  call timestamp(151,'done with Wannier functions')
endif
deallocate(norm)
!------------------------------!
! Hartree energy <W_n|V^H|W_n> !
!------------------------------!
allocate(ehart(nwann))
ehart=zzero
do n=1,nwann
  do ispn=1,nspinor
    ehart(n)=ehart(n)+lf_intgr_zdz(wanmt(1,1,1,1,ispn,n),wanir(1,1,ispn,n),&
      vhxcmt(1,1,1,1,ispn,n),vhxcir(1,1,ispn,n),(/0,0,0/),wanmt(1,1,1,1,ispn,n),&
      wanir(1,1,ispn,n))
  enddo
enddo
!-------------------------!
! XC potential and energy !
!-------------------------!
allocate(exc(nwann))
exc=zzero
call timer_start(12,reset=.true.)
call sic_genvxc(wanmt,wanir,vhxcmt,vhxcir,exc)
call timer_stop(12)
if (wproc) then
  write(151,'("time for XC potential : ",F8.3)')timer_get_value(12)
endif
sic_etot_correction=0.d0
ehart(:)=-1.d0*ehart(:)
exc(:)=-1.d0*exc(:)
do n=1,nwann
  sic_etot_correction=sic_etot_correction+dreal(ehart(n))+dreal(exc(n))
enddo
if (wproc) then
  write(151,*)
  write(151,'(2X,"wann",16X,"E_n^{H}",33X,"E_n^{XC}")')
  write(151,'(84("-"))')
  do n=1,nwann
    write(151,'(I4,4X,2G18.10,4X,2G18.10)')n,dreal(ehart(n)),dimag(ehart(n)),&
      dreal(exc(n)),dimag(exc(n))
  enddo
  write(151,*)
  write(151,'("Total energy correction : ",G18.10)')sic_etot_correction
  call flushifc(151)
endif
! multiply Wannier function by potential and change sign
do n=1,nwann
  do ispn=1,nspinor
    call lf_mult_zd(-zone,wanmt(1,1,1,1,ispn,n),wanir(1,1,ispn,n), &
      vhxcmt(1,1,1,1,ispn,n),vhxcir(1,1,ispn,n),wvmt(1,1,1,1,ispn,n),&
      wvir(1,1,ispn,n))    
  enddo
enddo
!----------------------------------!
! matrix elements of SIC potential !
!----------------------------------!
if (allocated(vwanme)) deallocate(vwanme)
allocate(vwanme(nmegqwan))
vwanme=zzero
! compute matrix elements of SIC potential
! vwan = <w_n|v_n|w_{n1,T}>
do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  v1l(:)=imegqwan(3:5,i)
  do ispn=1,nspinor    
    vwanme(i)=vwanme(i)+lf_dot_lf(.true.,wvmt(1,1,1,1,ispn,n),wvir(1,1,ispn,n),&
      v1l,wanmt(1,1,1,1,ispn,n1),wanir(1,1,ispn,n1))
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
      dreal(vwanme(i)),dimag(vwanme(i))
  enddo
endif
t1=0.d0
t2=0.d0
do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  v1l(:)=imegqwan(3:5,i)
  j=idxmegqwan(n1,n,-v1l(1),-v1l(2),-v1l(3))
  t1=t1+abs(vwanme(i)-dconjg(vwanme(j)))
  t2=max(t2,abs(vwanme(i)-dconjg(vwanme(j))))
enddo
if (wproc) then
  write(151,*)
  write(151,'("Maximum deviation from ""localization criterion"" : ",F12.6)')t2
  write(151,'("Average deviation from ""localization criterion"" : ",F12.6)')t1/nmegqwan
  write(151,*)
  write(151,'("Diagonal matrix elements")')
  write(151,'(2X,"wann",18X,"V_n")')
  write(151,'(44("-"))')
  do n=1,nwann
    j=idxmegqwan(n,n,0,0,0)
    write(151,'(I4,4X,2G18.10)')n,dreal(vwanme(j)),dimag(vwanme(j))
  enddo  
  call flushifc(151)
endif
if (wproc) then
  inquire(file="sic.hdf5",exist=exist)
  if (exist) then
    allocate(vwanme_old(nmegqwan))
    call hdf5_read("sic.hdf5","/","vwanme",vwanme_old(1),(/nmegqwan/))
    t1=0.d0
    do i=1,nmegqwan
      t1=t1+abs(vwanme(i)-vwanme_old(i))**2
    enddo
    t1=sqrt(t1/nmegqwan)
    write(151,*)
    write(151,'("SIC matrix elements RMS difference :",G18.10)')t1
    deallocate(vwanme_old)
  endif
endif
call sic_writevwan
if (wproc) close(151)
deallocate(vwanme)
deallocate(ehart)
deallocate(exc)
deallocate(vhxcmt,vhxcir)
deallocate(wanmt,wanir)
return
end




!subroutine test_lf(wmt,vmt,zsum1,zsum2)
!use modmain
!implicit none
!
!complex(8), intent(out) :: wmt(lmmaxvr,nrmtmax,natmtot)
!real(8), intent(out) :: vmt(lmmaxvr,nrmtmax,natmtot)
!complex(8), intent(out) :: zsum1
!complex(8), intent(out) :: zsum2
!complex(8), allocatable :: zfmt(:,:,:)
!complex(8), allocatable :: zgmt(:,:,:)
!real(8), allocatable :: dfmt(:,:,:)
!complex(8), allocatable :: zfmt_(:,:),zgmt_(:,:)
!real(8), allocatable :: dfmt_(:,:)
!integer ias,ir,lm,lm1,lm2,lm3
!real(8) d1,d2(2)
!complex(8) zt1,zt2
!complex(8) zf1(nrmtmax)
!complex(8), external :: gauntyry
!complex(8), external :: zfmtinp_
!
!
!allocate(zfmt(lmmaxvr,nrmtmax,natmtot))
!allocate(zgmt(lmmaxvr,nrmtmax,natmtot))
!allocate(dfmt(lmmaxvr,nrmtmax,natmtot))
!
!do ias=1,natmtot
!  do ir=1,nrmtmax
!    do lm=1,lmmaxvr
!      call random_number(d2)
!      zfmt(lm,ir,ias)=dcmplx(d2(1),d2(2))*10
!      call random_number(d1)
!      dfmt(lm,ir,ias)=d1*10
!    enddo
!  enddo
!enddo
!wmt=zfmt
!vmt=dfmt
!
!allocate(zfmt_(nrmtmax,lmmaxvr))
!allocate(zgmt_(nrmtmax,lmmaxvr))
!allocate(dfmt_(nrmtmax,lmmaxvr))
!
!! muffin-tin contribution
!do ias=1,natmtot
!  do lm1=1,lmmaxvr
!    zfmt_(:,lm1)=zfmt(lm1,:,ias)
!    dfmt_(:,lm1)=dfmt(lm1,:,ias)
!  enddo
!  do lm1=1,lmmaxvr
!    do lm2=1,lmmaxvr
!      do lm3=1,lmmaxvr
!        zt1=gauntyry(lm2l(lm1),lm2l(lm2),lm2l(lm3),&
!                 lm2m(lm1),lm2m(lm2),lm2m(lm3))
!        if (abs(zt1).gt.1d-12) then
!          do ir=1,nrmt(ias2is(ias))
!            zf1(ir)=dconjg(zfmt_(ir,lm1))*dfmt_(ir,lm2)*zfmt_(ir,lm3)*&
!              spr(ir,ias2is(ias))**2
!          enddo
!          zt2=zzero
!          do ir=1,nrmt(ias2is(ias))-1
!            zt2=zt2+0.5d0*(spr(ir+1,ias2is(ias))-spr(ir,ias2is(ias)))*&
!              (zf1(ir)+zf1(ir+1))
!          enddo
!          zsum1=zsum1+zt2*zt1
!        endif
!      enddo
!    enddo
!  enddo
!enddo !ias
!
!do ias=1,natmtot
!  do lm1=1,lmmaxvr
!    zfmt_(:,lm1)=zfmt(lm1,:,ias)
!    dfmt_(:,lm1)=dfmt(lm1,:,ias)
!    zgmt_(:,lm1)=zzero
!  enddo
!  do lm1=1,lmmaxvr
!    do lm2=1,lmmaxvr
!      do lm3=1,lmmaxvr
!        zt1=gauntyry(lm2l(lm3),lm2l(lm2),lm2l(lm1),&
!          lm2m(lm3),lm2m(lm2),lm2m(lm1))
!        write(180,*)zt1
!        if (abs(zt1).gt.1d-12) then
!          do ir=1,nrmt(ias2is(ias))
!            zgmt_(ir,lm3)=zgmt_(ir,lm3)+zfmt_(ir,lm1)*dfmt_(ir,lm2)*zt1
!          enddo
!        endif
!      enddo
!    enddo
!  enddo
!  !write(180,*)'ias=',ias
!  !write(180,*)zgmt_
!  do lm3=1,lmmaxvr
!    zgmt(lm3,:,ias)=zgmt_(:,lm3)
!  enddo
!  zsum2=zsum2+zfmtinp_(.true.,lmaxvr,nrmt(ias2is(ias)),spr(:,ias2is(ias)),&
!    lmmaxvr,zgmt(1,1,ias),zgmt(1,1,ias))
!enddo !ias
!
!!if (abs(zsum1-zsum2).gt.1d-8) 
!!write(*,*)abs(zsum1-zsum2),abs(zsum1-zsum2)/abs(zsum1),abs(zsum1-zsum2)/abs(zsum2)
!
!deallocate(zfmt,zgmt,dfmt,zfmt_,zgmt_,dfmt_)
!
!end


