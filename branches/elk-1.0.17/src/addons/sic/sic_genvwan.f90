subroutine sic_genvwan
use modmain
use mod_lf
use mod_nrkp
use modxcifc
use mod_addons_q
use mod_hdf5
implicit none
integer n,sz,i,j,n1,ispn,vtrl(3),nloc,n1loc,h,h1
real(8) t1,t2,vtrc(3)
integer vl(3)
! Wannier functions
complex(8), allocatable :: wanmt0(:,:,:,:,:)
complex(8), allocatable :: wanir0(:,:,:)
complex(8), allocatable :: vwanme_old(:)
complex(8), allocatable :: ene(:,:)
complex(8), allocatable :: vwank(:,:)
complex(8) z1
logical exist
integer n2,ik

! mpi grid layout
!          (2)
!     +----+----+--> T-vectos 
!     |    |    |
!     +----+----+--
! (1) |    |    |
!     +----+----+--
!     |    |    |
!     v
!  k-points and 
! Wannier functions

sic=.true.

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
call sic_wan(151)
allocate(ene(4,nwann))
call sic_pot(151,ene)
!----------------------------------!
! matrix elements of SIC potential !
!----------------------------------!
if (allocated(vwanme)) deallocate(vwanme)
allocate(vwanme(nmegqwan))
vwanme=zzero
allocate(wanmt0(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor))
allocate(wanir0(ngrtot,ntrloc,nspinor))
! compute matrix elements of SIC potential
! vwan = <w_n|v_n|w_{n1,T}>
do i=1,nmegqwan
  n=imegqwan(1,i)
  nloc=mpi_grid_map(nwann,dim_k,x=h,glob=n)
  n1=imegqwan(2,i)
  n1loc=mpi_grid_map(nwann,dim_k,x=h1,glob=n1)
  vl(:)=imegqwan(3:5,i)
  if (mpi_grid_x(dim_k).eq.h1) then
    call mpi_grid_send(wanmt(1,1,1,1,1,n1loc),lmmaxvr*nrmtmax*natmtot*ntrloc*nspinor,&
      dims=(/dim_k/),dest=(/h/),tag=1)
    call mpi_grid_send(wanir(1,1,1,n1loc),ngrtot*ntrloc*nspinor,&
      dims=(/dim_k/),dest=(/h/),tag=2)
  endif
  if (mpi_grid_x(dim_k).eq.h) then
    call mpi_grid_recieve(wanmt0(1,1,1,1,1),lmmaxvr*nrmtmax*natmtot*ntrloc*nspinor,&
      dims=(/dim_k/),src=(/h1/),tag=1)
    call mpi_grid_recieve(wanir0(1,1,1),ngrtot*ntrloc*nspinor,&
      dims=(/dim_k/),src=(/h1/),tag=2)
  endif
  if (mpi_grid_x(dim_k).eq.h) then
    do ispn=1,nspinor    
      vwanme(i)=vwanme(i)+lf_dot_lf(.true.,wvmt(1,1,1,1,ispn,nloc),wvir(1,1,ispn,nloc),&
        vl,wanmt0(1,1,1,1,ispn),wanir0(1,1,ispn))
    enddo
  endif
enddo
call mpi_grid_reduce(vwanme(1),nmegqwan,dims=(/dim_k/))
t1=0.d0
t2=0.d0
do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  vl(:)=imegqwan(3:5,i)
  j=idxmegqwan(n1,n,-vl(1),-vl(2),-vl(3))
  t1=t1+abs(vwanme(i)-dconjg(vwanme(j)))
  t2=max(t2,abs(vwanme(i)-dconjg(vwanme(j))))
enddo
if (wproc) then
  call timestamp(151,"done with matrix elements")
  write(151,*)
  write(151,'("Number of Wannier transitions : ",I6)')nmegqwan
  write(151,'("Matrix elements of SIC potential &
    &(n n1  T  <w_n|v_n|w_{n1,T}>)")')
  do i=1,nmegqwan
    write(151,'(I4,4X,I4,4X,3I3,4X,2G18.10)')imegqwan(:,i),&
      dreal(vwanme(i)),dimag(vwanme(i))
  enddo
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
if (wproc) write(151,*)
! check hermiticity of V_nn'(k)
allocate(vwank(nwann,nwann))
vwank=zzero
do ik=1,nkpt
  do i=1,nmegqwan
    n1=imegqwan(1,i)
    n2=imegqwan(2,i)
    vtrl(:)=imegqwan(3:5,i)
    vtrc(:)=vtrl(1)*avec(:,1)+vtrl(2)*avec(:,2)+vtrl(3)*avec(:,3)
    z1=exp(zi*dot_product(vkc(:,ik),vtrc(:)))
    vwank(n1,n2)=vwank(n1,n2)+z1*vwanme(i)
  enddo
  t1=0.d0
  do n1=1,nwann
    do n2=1,nwann
      t1=max(t1,abs(vwank(n1,n2)-dconjg(vwank(n2,n1))))
    enddo
  enddo
  if (wproc) then
    write(151,'("ik : ",I4,"   max.herm.err : ",G18.10 )')ik,t1
  endif
enddo
deallocate(vwank)
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
deallocate(ene)
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


