subroutine sic_main
use modmain
use mod_nrkp
use mod_hdf5
use mod_sic
implicit none
integer n,sz,i,j,i1,j1,n1,ispn,vtrl(3)
real(8) t1,t2,t3,vtrc(3)
integer vl(3)
! Wannier functions
complex(8), allocatable :: vwanme_old(:)
complex(8), allocatable :: ene(:,:)
complex(8), allocatable :: vwank(:,:)
complex(8) z1
logical exist
integer n2,ik

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
  sz=lmmaxvr*nmtloc+ngrloc
  sz=16*sz*ntr*nspinor*(2*nwann+2)/1024/1024
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
call sic_wan(151)
allocate(ene(4,nwann))
call sic_pot(151,ene)
!----------------------------------!
! matrix elements of SIC potential !
!----------------------------------!
allocate(vwanme_old(nmegqwan))
vwanme_old=vwanme
vwanme=zzero
! compute matrix elements of SIC potential
!  vwanme = <w_n|v_n|w_{n1,T}>
do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  vl(:)=imegqwan(3:5,i)
  do ispn=1,nspinor    
    vwanme(i)=vwanme(i)+sic_dot_ll(wvmt(1,1,1,ispn,n),wvir(1,1,ispn,n),&
      wanmt(1,1,1,ispn,n1),wanir(1,1,ispn,n1),vl,twanmt(1,1,n),twanmt(1,1,n1))
  enddo
enddo
t1=0.d0
t2=-1.d0
t3=0.d0
do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  vl(:)=imegqwan(3:5,i)
  j=idxmegqwan(n1,n,-vl(1),-vl(2),-vl(3))
  t1=t1+abs(vwanme(i)-dconjg(vwanme(j)))
  if (abs(vwanme(i)-dconjg(vwanme(j))).ge.t2) then
    t2=abs(vwanme(i)-dconjg(vwanme(j)))
    i1=i
    j1=j
  endif
  t3=t3+abs(vwanme(i)-vwanme_old(i))**2
enddo
if (wproc) then
  call timestamp(151,"done with matrix elements")
  write(151,*)
  write(151,'("Number of Wannier transitions : ",I6)')nmegqwan
  write(151,'("Matrix elements of SIC potential &
    &(n n1  <w_n|v_n|w_n1}>)")')
  do i=1,nmegqwan
    vl(:)=imegqwan(3:5,i)
    if (all(vl.eq.0)) then
      write(151,'(I4,4X,I4,4X,2G18.10)')imegqwan(1:2,i),&
        dreal(vwanme(i)),dimag(vwanme(i))
    endif
  enddo
  write(151,*)
  write(151,'("Maximum deviation from ""localization criterion"" : ",F12.6)')t2
  write(151,'("Average deviation from ""localization criterion"" : ",F12.6)')t1/nmegqwan
  write(151,*)
  write(151,'("Matrix elements with maximum difference : ",2I6)')i1,j1
  write(151,'(I4,4X,I4,4X,3I4,4X,2G18.10)')imegqwan(:,i1),&
        dreal(vwanme(i1)),dimag(vwanme(i1))
  write(151,'(I4,4X,I4,4X,3I4,4X,2G18.10)')imegqwan(:,j1),&
        dreal(vwanme(j1)),dimag(vwanme(j1))
  write(151,*)
  write(151,'("Diagonal matrix elements")')
  write(151,'(2X,"wann",18X,"V_n")')
  write(151,'(44("-"))')
  do n=1,nwann
    j=idxmegqwan(n,n,0,0,0)
    write(151,'(I4,4X,2G18.10)')n,dreal(vwanme(j)),dimag(vwanme(j))
  enddo  
  t3=sqrt(t3/nmegqwan)
  write(151,*)
  write(151,'("SIC matrix elements RMS difference :",G18.10)')t3  
  call flushifc(151)
endif
deallocate(vwanme_old)
! check hermiticity of V_nn'(k)
allocate(vwank(nwann,nwann))
do ik=1,nkpt
  vwank=zzero
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
    write(151,*)
    write(151,'("ik : ",I4,"   max.herm.err : ",G18.10 )')ik,t1
!    write(151,*)
!    do n1=1,nwann
!      write(151,'(5X,255F12.7)')(dreal(vwank(n1,n2)),n2=1,nwann)
!    enddo
!    write(151,*)
!    do n1=1,nwann
!      write(151,'(5X,255F12.7)')(dimag(vwank(n1,n2)),n2=1,nwann)
!    enddo
  endif
enddo
deallocate(vwank)
!if (wproc) then
!  inquire(file="sic.hdf5",exist=exist)
!  if (exist) then
!    allocate(vwanme_old(nmegqwan))
!    call hdf5_read("sic.hdf5","/","vwanme",vwanme_old(1),(/nmegqwan/))
!    t1=0.d0
!    do i=1,nmegqwan
!      t1=t1+abs(vwanme(i)-vwanme_old(i))**2
!    enddo
!    t1=sqrt(t1/nmegqwan)
!    write(151,*)
!    write(151,'("SIC matrix elements RMS difference :",G18.10)')t1
!    deallocate(vwanme_old)
!  endif
!endif
if (wproc) close(151)
! flag that now we have computed sic potential and wannier functions
tsic_wv=.true.
! write to HDF5 file after last iteration
if (isclsic.eq.nsclsic) call sic_writevwan
deallocate(ene)
return
end
!
!
!
!
!!subroutine test_lf(wmt,vmt,zsum1,zsum2)
!!use modmain
!!implicit none
!!
!!complex(8), intent(out) :: wmt(lmmaxvr,nrmtmax,natmtot)
!!real(8), intent(out) :: vmt(lmmaxvr,nrmtmax,natmtot)
!!complex(8), intent(out) :: zsum1
!!complex(8), intent(out) :: zsum2
!!complex(8), allocatable :: zfmt(:,:,:)
!!complex(8), allocatable :: zgmt(:,:,:)
!!real(8), allocatable :: dfmt(:,:,:)
!!complex(8), allocatable :: zfmt_(:,:),zgmt_(:,:)
!!real(8), allocatable :: dfmt_(:,:)
!!integer ias,ir,lm,lm1,lm2,lm3
!!real(8) d1,d2(2)
!!complex(8) zt1,zt2
!!complex(8) zf1(nrmtmax)
!!complex(8), external :: gauntyry
!!complex(8), external :: zfmtinp_
!!
!!
!!allocate(zfmt(lmmaxvr,nrmtmax,natmtot))
!!allocate(zgmt(lmmaxvr,nrmtmax,natmtot))
!!allocate(dfmt(lmmaxvr,nrmtmax,natmtot))
!!
!!do ias=1,natmtot
!!  do ir=1,nrmtmax
!!    do lm=1,lmmaxvr
!!      call random_number(d2)
!!      zfmt(lm,ir,ias)=dcmplx(d2(1),d2(2))*10
!!      call random_number(d1)
!!      dfmt(lm,ir,ias)=d1*10
!!    enddo
!!  enddo
!!enddo
!!wmt=zfmt
!!vmt=dfmt
!!
!!allocate(zfmt_(nrmtmax,lmmaxvr))
!!allocate(zgmt_(nrmtmax,lmmaxvr))
!!allocate(dfmt_(nrmtmax,lmmaxvr))
!!
!!! muffin-tin contribution
!!do ias=1,natmtot
!!  do lm1=1,lmmaxvr
!!    zfmt_(:,lm1)=zfmt(lm1,:,ias)
!!    dfmt_(:,lm1)=dfmt(lm1,:,ias)
!!  enddo
!!  do lm1=1,lmmaxvr
!!    do lm2=1,lmmaxvr
!!      do lm3=1,lmmaxvr
!!        zt1=gauntyry(lm2l(lm1),lm2l(lm2),lm2l(lm3),&
!!                 lm2m(lm1),lm2m(lm2),lm2m(lm3))
!!        if (abs(zt1).gt.1d-12) then
!!          do ir=1,nrmt(ias2is(ias))
!!            zf1(ir)=dconjg(zfmt_(ir,lm1))*dfmt_(ir,lm2)*zfmt_(ir,lm3)*&
!!              spr(ir,ias2is(ias))**2
!!          enddo
!!          zt2=zzero
!!          do ir=1,nrmt(ias2is(ias))-1
!!            zt2=zt2+0.5d0*(spr(ir+1,ias2is(ias))-spr(ir,ias2is(ias)))*&
!!              (zf1(ir)+zf1(ir+1))
!!          enddo
!!          zsum1=zsum1+zt2*zt1
!!        endif
!!      enddo
!!    enddo
!!  enddo
!!enddo !ias
!!
!!do ias=1,natmtot
!!  do lm1=1,lmmaxvr
!!    zfmt_(:,lm1)=zfmt(lm1,:,ias)
!!    dfmt_(:,lm1)=dfmt(lm1,:,ias)
!!    zgmt_(:,lm1)=zzero
!!  enddo
!!  do lm1=1,lmmaxvr
!!    do lm2=1,lmmaxvr
!!      do lm3=1,lmmaxvr
!!        zt1=gauntyry(lm2l(lm3),lm2l(lm2),lm2l(lm1),&
!!          lm2m(lm3),lm2m(lm2),lm2m(lm1))
!!        write(180,*)zt1
!!        if (abs(zt1).gt.1d-12) then
!!          do ir=1,nrmt(ias2is(ias))
!!            zgmt_(ir,lm3)=zgmt_(ir,lm3)+zfmt_(ir,lm1)*dfmt_(ir,lm2)*zt1
!!          enddo
!!        endif
!!      enddo
!!    enddo
!!  enddo
!!  !write(180,*)'ias=',ias
!!  !write(180,*)zgmt_
!!  do lm3=1,lmmaxvr
!!    zgmt(lm3,:,ias)=zgmt_(:,lm3)
!!  enddo
!!  zsum2=zsum2+zfmtinp_(.true.,lmaxvr,nrmt(ias2is(ias)),spr(:,ias2is(ias)),&
!!    lmmaxvr,zgmt(1,1,ias),zgmt(1,1,ias))
!!enddo !ias
!!
!!!if (abs(zsum1-zsum2).gt.1d-8) 
!!!write(*,*)abs(zsum1-zsum2),abs(zsum1-zsum2)/abs(zsum1),abs(zsum1-zsum2)/abs(zsum2)
!!
!!deallocate(zfmt,zgmt,dfmt,zfmt_,zgmt_,dfmt_)
!!
!!end
!
!
