subroutine sic_genvwan
use modmain
use mod_lf
use mod_nrkp
use modxcifc
use mod_addons_q
use mod_hdf5
implicit none

integer n
integer ik,ikloc,j,sz,i,itp
integer n1,n2,ispn,irloc,nrmtloc
integer itr,it,itloc,ir,m,ias
real(8) t1,t2
integer v1l(3),lm1,lm2,lm3
character*12 c1,c2,c3
character*100 path
integer ntp
real(8), allocatable :: tp(:,:)
complex(8), allocatable :: ylm(:,:)
complex(8) zt1
real(8), external :: gaunt

! arrays for Wannier functions
complex(8), allocatable :: wanmt0(:,:,:,:,:)
complex(8), allocatable :: wanir0(:,:,:)

real(8), allocatable :: rhowanir(:)
complex(8), allocatable :: f1mt(:,:)
complex(8), allocatable :: f2mt(:,:)
complex(8), allocatable :: f3mt(:,:)
complex(8), allocatable :: f4mt(:,:,:)

complex(8), allocatable :: zm1(:,:,:)
real(8), allocatable :: f3(:),f4(:),f5(:)

complex(8), allocatable :: wfmt(:,:,:,:)
complex(8), allocatable :: wfir(:,:)
complex(8), allocatable :: a2(:,:)
real(8) v2(3),v3(3)
complex(8) expikr
integer i1,i2,i3,is,io,ig

integer lm
real(8), allocatable :: vx(:),vc(:)
complex(8) z1,expikt
real(8) vtrc(3)

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
! read Fermi energy
if (mpi_grid_root()) call readfermi
call mpi_grid_bcast(efermi)

call lf_init(lf_maxt,dim2)
wproc=mpi_grid_root()
if (wproc) then
  open(151,file='SIC.OUT',form='FORMATTED',status='REPLACE')
endif
call genwfnr(151,.false.)  
if (wproc) then
  call timestamp(151,'done with wavefunctions')
  call flushifc(151)
endif
! get all Wannier transitions
all_wan_ibt=.true.
call getimegqwan(all_wan_ibt)
! compute Fourier transform of <n,T=0|H^{LDA}|n',T'>
allocate(zm1(nwann,nwann,nkptnrloc))
zm1=zzero
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  do n1=1,nwann
    do n2=1,nwann
      do j=1,nstsv
        zm1(n1,n2,ikloc)=zm1(n1,n2,ikloc)+dconjg(wann_c(n1,j,ikloc))*&
          wann_c(n2,j,ikloc)*evalsvnr(j,ik)
      enddo
    enddo
  enddo
enddo 
! compute <n,T=0|H^{LDA}|n',T'>
if (allocated(hwan)) deallocate(hwan)
allocate(hwan(nmegqwan))
hwan=zzero
do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  v1l(:)=imegqwan(3:5,i)
  vtrc(:)=v1l(1)*avec(:,1)+v1l(2)*avec(:,2)+v1l(3)*avec(:,3)
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    expikt=exp(-1.d0*zi*dot_product(vkcnr(:,ik),vtrc(:)))
    hwan(i)=hwan(i)+expikt*zm1(n,n1,ikloc)
  enddo
enddo
call mpi_grid_reduce(hwan(1),nmegqwan,dims=(/dim_k/))
hwan(:)=hwan(:)/nkptnr
deallocate(zm1)


if (allocated(vwanmt)) deallocate(vwanmt)
allocate(vwanmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor,nwann))
if (allocated(vwanir)) deallocate(vwanir)
allocate(vwanir(ngrtot,ntrloc,nspinor,nwann))
vwanmt=zzero
vwanir=zzero

call sic_genvhart

! restore wproc
wproc=mpi_grid_root()
if (wproc) then
  write(151,*)
  write(151,'("time for q-vectors : ",F8.3)')timer_get_value(10)
  write(151,'("time for Hartree potential : ",F8.3)')timer_get_value(11)
endif


! generate Wannier functions on a mesh
if (allocated(wanmt)) deallocate(wanmt)
allocate(wanmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor,nwann))
if (allocated(wanir)) deallocate(wanir)
allocate(wanir(ngrtot,ntrloc,nspinor,nwann))
if (wproc) then
  sz=lmmaxvr*nrmtmax*natmtot+ngrtot
  sz=16*sz*nspinor*nwann*ntrloc/1024/1024
  write(151,*)
  write(151,'("Size of real-space Wannier functions arrays (MB) : ",I6)')sz
  write(151,*)
  call flushifc(151)
endif
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
  write(151,'("MT part : ",F8.3)')timer_get_value(1)
  write(151,'("IT part : ",F8.3)')timer_get_value(2)
  call flushifc(151)
endif
if (wproc) then
  call timestamp(151,'done with Wannier functions')
endif


! test
allocate(wfmt(lmmaxvr,nrmtmax,natmtot,nspinor))
allocate(wfir(ngrtot,nspinor))
allocate(a2(nwann,nstsv))
t1=0.d0
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  a2=zzero
  do j=1,nstsv
    wfmt=zzero
    wfir=zzero
    do ispn=1,nspinor
      do ias=1,natmtot
        is=ias2is(ias)
        do ir=1,nrmt(is)
          do lm=1,lmmaxvr
            do io=1,nufr(lm2l(lm),is)
              wfmt(lm,ir,ias,ispn)=wfmt(lm,ir,ias,ispn)+&
                ufr(ir,lm2l(lm),io,ias)*wfsvmtloc(lm,io,ias,ispn,j,ikloc)
            enddo
          enddo
        enddo
      enddo !ias
      do ig=1,ngknr(ikloc)
        wfir(igfft(igkignr(ig,ikloc)),ispn)=wfsvitloc(ig,ispn,j,ikloc)
      enddo
      call zfftifc(3,ngrid,1,wfir(:,ispn))
      wfir(:,ispn)=wfir(:,ispn)/sqrt(omega)
    enddo !ispn
    ir=0
    do i3=0,ngrid(3)-1
      v2(3)=dble(i3)/dble(ngrid(3))
      do i2=0,ngrid(2)-1
        v2(2)=dble(i2)/dble(ngrid(2))
        do i1=0,ngrid(1)-1
          v2(1)=dble(i1)/dble(ngrid(1))
          ir=ir+1
          call r3mv(avec,v2,v3)
          expikr=exp(zi*dot_product(vkcnr(:,ik),v3(:)))
          wfir(ir,:)=expikr*wfir(ir,:)
        enddo
      enddo
    enddo
    do n=1,nwann
      do ispn=1,nspinor
        a2(n,j)=a2(n,j)+lf_dotblh(.true.,vkcnr(:,ik),wanmt(1,1,1,1,ispn,n),&
          wanir(1,1,ispn,n),wfmt(1,1,1,ispn),wfir(1,ispn))
      enddo
    enddo
  enddo !j
  do n=1,nwann
    do j=1,nstsv
      t1=max(t1,abs(wann_c(n,j,ikloc)-dconjg(a2(n,j))))
    enddo
  enddo
enddo !ikloc
deallocate(wfmt,wfir,a2)
call mpi_grid_reduce(t1,dims=(/dim_k/),op=op_max)
if (wproc) then
  write(151,*)
  write(151,'("Maximum deviation from exact expansion : ",G18.10)')t1
endif
! deallocate unnecessary arrays
deallocate(wfsvmtloc)
deallocate(wfsvitloc)
deallocate(evecfvloc)
deallocate(evecsvloc)
deallocate(wann_c)


call timer_start(12,reset=.true.)
ntp=1000
allocate(tp(2,ntp))
allocate(ylm(lmmaxvr,ntp))
call sphcover(ntp,tp)
do itp=1,ntp 
  call genylm(lmaxvr,tp(1,itp),ylm(1,itp))
enddo
m=max(ntp,ngrtot)
allocate(rhowanir(ngrtot))
allocate(f3(m),f4(m),f5(m))
allocate(vx(m),vc(m))
allocate(f4mt(lmmaxvr,nrmtmax,natmtot))
! add XC potential to Coulomb
do n=1,nwann
  do itloc=1,ntrloc
    f4mt=zzero
! muffin-tin part
    do ias=1,natmtot
      nrmtloc=mpi_grid_map(nrmt(ias2is(ias)),dim_k)
      !do ir=1,nrmt(ias2is(ias))
      do irloc=1,nrmtloc
        ir=mpi_grid_map(nrmt(ias2is(ias)),dim_k,loc=irloc)
! compute charge density on a sphere
        f5=0.d0
        do itp=1,ntp
          do ispn=1,nspinor
            zt1=zzero            
            do lm=1,lmmaxvr
              zt1=zt1+wanmt(lm,ir,ias,itloc,ispn,n)*ylm(lm,itp)
            enddo
            f5(itp)=f5(itp)+abs(zt1)**2
          enddo
        enddo !itp
        call xcifc(xctype,n=ntp,rho=f5,ex=f3,ec=f4,vx=vx,vc=vc)
! save XC potential
        f5(1:ntp)=vx(1:ntp)+vc(1:ntp)
! expand XC potential in spherical harmonics
        do lm=1,lmmaxvr
          zt1=zzero
          do itp=1,ntp
            zt1=zt1+dconjg(ylm(lm,itp))*f5(itp)
          enddo
          f4mt(lm,ir,ias)=fourpi*zt1/ntp
        enddo !lm
      enddo !irloc
    enddo !ias
    call mpi_grid_reduce(f4mt(1,1,1),lmmaxvr*nrmtmax*natmtot,dims=(/dim_k/),&
      all=.true.)
    vwanmt(:,:,:,itloc,1,n)=vwanmt(:,:,:,itloc,1,n)+f4mt(:,:,:)
    rhowanir(:)=dreal(dconjg(wanir(:,itloc,1,n))*wanir(:,itloc,1,n))
    if (spinpol) then
      rhowanir(:)=rhowanir(:)+&
        dreal(dconjg(wanir(:,itloc,2,n))*wanir(:,itloc,2,n))
    endif
    call xcifc(xctype,n=ngrtot,rho=rhowanir(:),ex=f3,ec=f4,vx=vx,vc=vc)
    vwanir(:,itloc,1,n)=vwanir(:,itloc,1,n)+vc(1:ngrtot)+vx(1:ngrtot)
  enddo
enddo
deallocate(vx,vc,f3,f4,rhowanir,f4mt)
call timer_stop(12)
if (wproc) then
  write(151,'("time for XC potential : ",F8.3)')timer_get_value(12)
endif

if (spinpol) then
  do n=1,nwann
    vwanmt(:,:,:,:,2,n)=vwanmt(:,:,:,:,1,n)
    vwanir(:,:,2,n)=vwanir(:,:,1,n)
  enddo
endif
call timer_start(13,reset=.true.)
! multiply potential by Wannier function and change sign
allocate(f1mt(nrmtmax,lmmaxvr))
allocate(f2mt(nrmtmax,lmmaxvr))
allocate(f3mt(nrmtmax,lmmaxvr))
do n=1,nwann
  do itloc=1,ntrloc
    do ispn=1,nspinor     
      do ias=1,natmtot
        f3mt=zzero
        do lm1=1,lmmaxvr
          f1mt(:,lm1)=vwanmt(lm1,:,ias,itloc,ispn,n)
          f2mt(:,lm1)=wanmt(lm1,:,ias,itloc,ispn,n)
        enddo
        do lm1=1,lmmaxvr
          do lm2=1,lmmaxvr
            do lm3=1,lmmaxvr
              t1=gaunt(lm2l(lm3),lm2l(lm1),lm2l(lm2),&
                       lm2m(lm3),lm2m(lm1),lm2m(lm2))
              if (abs(t1).gt.1d-8) then
                do ir=1,nrmt(ias2is(ias))
                  f3mt(ir,lm3)=f3mt(ir,lm3)+f1mt(ir,lm1)*f2mt(ir,lm2)*t1
                enddo
              endif
            enddo
          enddo
        enddo
        do lm3=1,lmmaxvr
          vwanmt(lm3,:,ias,itloc,ispn,n)=-f3mt(:,lm3)
        enddo
      enddo !ias
      vwanir(:,itloc,ispn,n)=-vwanir(:,itloc,ispn,n)*wanir(:,itloc,ispn,n)
    enddo !ispn
  enddo !itloc
enddo !n
deallocate(f1mt,f2mt,f3mt)
call timer_stop(13)
if (wproc) then
  write(151,'("time for V*WF product : ",F8.3)')timer_get_value(13)
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
  write(151,'("   n    Re H_nn     Im H_nn     Re V_n      Im V_n")')
  write(151,'(70("-"))')
  do n=1,nwann
    j=idxmegqwan(n,n,0,0,0)
    write(151,'(I4,4F12.6)')n,dreal(hwan(j)),dimag(hwan(j)),&
      dreal(vwan(j)),dimag(vwan(j))
  enddo  
  call flushifc(151)
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
  call hdf5_write("sic.hdf5","/","hwan",hwan(1),(/nmegqwan/))
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
deallocate(hwan,vwan)
deallocate(wanmt,wanir,vwanmt,vwanir)
return
end

