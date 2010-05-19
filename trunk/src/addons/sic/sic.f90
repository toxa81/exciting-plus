subroutine sic
use modmain
use mod_lfa
use mod_nrkp
use modxcifc
implicit none

integer i1,i2,i3,n
integer ik,ikloc,ik1,j,ig,sz,i,isym
integer n1,n2,ispn
integer itr,ntrloc,it,jt,itloc,ir,m,ias
real(8) vtrc(3),t1,t2
integer v1l(3)
complex(8) expikt

! arrays for Wannier functions
complex(8), allocatable :: wanmt(:,:,:,:,:,:)
complex(8), allocatable :: wanir(:,:,:,:)
complex(8), allocatable :: wanmt0(:,:,:,:,:)
complex(8), allocatable :: wanir0(:,:,:)
! Hartree potential of a Wanier state
complex(8), allocatable :: vhwanmt(:,:,:,:,:)
complex(8), allocatable :: vhwanir(:,:,:)

complex(8), allocatable :: megqwan1(:,:,:) 

complex(8), allocatable :: rhowanmt(:,:,:,:,:)
complex(8), allocatable :: rhowanir(:,:,:)
complex(8), allocatable :: identmt(:,:,:,:)
complex(8), allocatable :: identir(:,:)
complex(8), allocatable :: f1mt(:,:,:,:)
complex(8), allocatable :: f1ir(:,:)
complex(8), allocatable :: f2mt(:,:,:)
complex(8), allocatable :: f2ir(:)
complex(8), allocatable :: rhokwanmt(:,:,:)
complex(8), allocatable :: rhokwanir(:)
complex(8), allocatable :: vckwanmt(:,:,:)
complex(8), allocatable :: vckwanir(:)
complex(8), allocatable :: vcwanmt(:,:,:,:,:)
complex(8), allocatable :: vcwanir(:,:,:)
real(8) spzn1(maxspecies)
complex(8), allocatable :: vsic(:,:,:)
complex(8), allocatable :: ubare(:,:,:)
complex(8), allocatable :: h0wan(:,:,:)
complex(8), allocatable :: zm1(:,:,:)
real(8), allocatable :: f3(:),f4(:)

integer nvq0loc,iqloc,iq
real(8), allocatable :: vx(:),vc(:)
integer idm
complex(8), allocatable :: ovlm(:,:,:)
complex(8), external :: zfinp_
complex(8) z1,z2
integer np
real(8), allocatable :: vpl(:,:)
real(8), allocatable :: fp(:)
character*100 qnm


! mpi grid layout
!          (3)
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
wproc=mpi_grid_root()
! read the density and potentials from file
call readstate
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
call geturf
call genurfprod
! read Fermi energy
if (mpi_grid_root()) call readfermi
call mpi_grid_bcast(efermi)

if (wproc) then
  open(151,file='SIC.OUT',form='FORMATTED',status='REPLACE')
endif

call lfa_init(1)
call genwfnr(151,.true.)  
call init_qbz(.true.)
if (spinpol) then
  if (allocated(spinor_ud)) deallocate(spinor_ud)
  allocate(spinor_ud(2,nstsv,nkptnr))
  spinor_ud=0
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do j=1,nstsv
      t1=sum(abs(evecsvloc(1:nstfv,j,ikloc)))
      if (t1.gt.1d-10) spinor_ud(1,j,ik)=1
      t1=sum(abs(evecsvloc(nstfv+1:nstsv,j,ikloc)))
      if (t1.gt.1d-10) spinor_ud(2,j,ik)=1
    enddo
  enddo
  call mpi_grid_reduce(spinor_ud(1,1,1),2*nstsv*nkptnr,dims=(/dim_k/),all=.true.)
endif  
! get all Wannier transitions
call getimegqwan(.true.)
if (wproc) then
  call timestamp(151,'done with wavefunctions')
endif

! create q-directories
if (mpi_grid_root()) then
  call system("mkdir -p q")
  do iq=1,nvq0
    call getqdir(iq,ivq0m_list(:,iq),qnm)
    call system("mkdir -p "//trim(qnm))
  enddo
endif
wannier_megq=.true.
all_wan_ibt=.true.
! distribute q-vectors along 3-rd dimention
nvq0loc=mpi_grid_map(nvq0,dim_q)
! distribute translations along 3-rd dimention
ntrloc=mpi_grid_map(ntr,dim_t)

! main loop over q-points
do iqloc=1,nvq0loc
  iq=mpi_grid_map(nvq0,dim_q,loc=iqloc)
  write(*,*)'iq=',iq
  call genmegq(iq,.true.)
! TODO: ngvecme is not known before genmegq. This way of allocation is
!  not good looking
  if (.not.allocated(megqwan1)) then
    allocate(megqwan1(nwann,ngvecme,nvq0))
    megqwan1=zzero
  endif
  do n=1,nwann
    megqwan1(n,:,iq)=megqwan(idxmegqwan(n,n,0,0,0),:)
  enddo
enddo
wproc=mpi_grid_root()
call mpi_grid_reduce(megqwan1(1,1,1),nwann*ngvecme*nvq0,dims=(/dim_q/), &
  all=.true.)
  
write(100+iproc,*)'megqwan1=',megqwan1

!allocate(h0wan(nwann,nwann,ntr))
!! compute <n,T=0|H^{LDA}|n',T'>
!allocate(zm1(nwann,nwann,nkptnrloc))
!zm1=zzero
!do ikloc=1,nkptnrloc
!  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
!  do n1=1,nwann
!    do n2=1,nwann
!      do j=1,nstsv
!        zm1(n1,n2,ikloc)=zm1(n1,n2,ikloc)+dconjg(wann_c(n1,j,ikloc))*&
!          wann_c(n2,j,ikloc)*evalsvnr(j,ik)
!      enddo
!    enddo
!  enddo
!enddo 
!h0wan=zzero
!do it=1,ntr
!  vtrc(:)=vtl(1,it)*avec(:,1)+vtl(2,it)*avec(:,2)+vtl(3,it)*avec(:,3)
!  do ikloc=1,nkptnrloc
!    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
!    expikt=exp(-1.d0*zi*dot_product(vkcnr(:,ik),vtrc(:)))
!    h0wan(:,:,it)=h0wan(:,:,it)+expikt*zm1(:,:,ikloc)
!  enddo
!enddo
!call mpi_grid_reduce(h0wan(1,1,1),nwann*nwann*ntr,dims=(/dim_k/))
!h0wan(:,:,:)=h0wan(:,:,:)/nkptnr
!!if (wproc) then
!!  write(151,*)
!!  write(151,'("hopping parameters <w_n1|H^{LDA}|w_{n2,T}>")')
!!  do it=1,ntr
!!    write(151,'("  translation : ",3I4)')vtl(:,it)
!!    write(151,'("    real part : ")')  
!!    do n1=1,nwann
!!      write(151,'(4X,255F12.6)')(dreal(h0wan(n1,n2,it)),n2=1,nwann)
!!    enddo
!!    write(151,'("    image part : ")')  
!!    do n1=1,nwann
!!      write(151,'(4X,255F12.6)')(dimag(h0wan(n1,n2,it)),n2=1,nwann)
!!    enddo
!!  enddo
!!endif
!if (wproc) then
!  call timestamp(151,'done with hoppings')
!endif
!! deallocate unnecessary wave-functions
!deallocate(wfsvmtloc)
!deallocate(wfsvitloc)
!deallocate(evecfvloc)
!deallocate(evecsvloc)
!deallocate(wann_c)

!if (wproc) then
!  write(151,*)
!  write(151,'("Number of translations : ",I4)')ntr
!endif

allocate(wanmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor,nwann))
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
  call gen_wann_func(vtl(1,itr),ngknr,vgkcnr,igkignr,wanmt0,wanir0)
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

if (wproc) then
  write(151,*)
  write(151,'("Number of Wannier transitions : ",I6)')nmegqwan
  write(151,'("List of Wannier transitions (n n1 T) ")')
  do i=1,nmegqwan
    write(151,'(I4,4X,I4,4X,3I3)')imegqwan(:,i)
  enddo
endif

! check orthonormality
t2=0.d0
do i=1,nmegqwan
  n=imegqwan(1,nmegqwan)
  n1=imegqwan(2,nmegqwan)
  v1l(:)=imegqwan(3:5,nmegqwan)
  z1=0
  do ispn=1,nspinor
    z1=z1+lfa_dotp(.true.,v1l,wanmt(1,1,1,1,ispn,n),wanir(1,1,ispn,n),&
      wanmt(1,1,1,1,ispn,n1),wanir(1,1,ispn,n1))
  enddo
  if (n.eq.n1.and.v1l(1).eq.0.and.v1l(2).eq.0.and.v1l(3).eq.0) then
    z1=z1-zone
  endif
  t2=max(t2,abs(z1))
enddo
if (wproc) then
  write(151,*)
  write(151,'("Maximum deviation from norm : ",F12.6)')t2
  call flushifc(151)
endif
call pstop


!allocate(vcwanmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nwann))
!allocate(vcwanir(ngrtot,ntrloc,nwann))
!
!vcwanmt=zzero
!vcwanir=zzero
!
!allocate(f2mt(lmmaxvr,nrmtmax,natmtot))
!allocate(f2ir(ngrtot))
!do ikloc=1,nkptnrloc
!  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
!  do n=1,nwann
!    rhokwanmt=zzero
!    rhokwanir=zzero
!! rho_{nk}=\sum_{T}exp^{+ikT}rho_n(r-T)=\sum_{T'}exp^{-ikT'}rho_n(r+T')    
!    do itloc=1,ntrloc
!      it=mpi_grid_map(ntr,dim_t,loc=itloc)
!      vtrc(:)=vtl(1,it)*avec(:,1)+vtl(2,it)*avec(:,2)+vtl(3,it)*avec(:,3)
!      expikt=exp(-1.d0*zi*dot_product(vkcnr(:,ik),vtrc(:)))
!      do ispn=1,nspinor
!        call lfa_sht('B',wanmt(1,1,1,itloc,ispn,n),f2mt)
!        rhokwanmt(:,:,:)=rhokwanmt(:,:,:)+expikt*dconjg(f2mt(:,:,:))*f2mt(:,:,:)
!        rhokwanir(:)=rhokwanir(:)+expikt*dconjg(wanir(:,itloc,ispn,n))*wanir(:,itloc,ispn,n)
!      enddo
!    enddo
!    call mpi_grid_reduce(rhokwanmt(1,1,1),lmmaxvr*nrmtmax*natmtot,dims=(/dim_t/),all=.true.)
!    call mpi_grid_reduce(rhokwanir(1),ngrtot,dims=(/dim_t/),all=.true.)
!    call lfa_sht('F',rhokwanmt,rhokwanmt)    
!    call potcoul_(rhokwanmt,rhokwanir,vckwanmt,vckwanir)
!! Vc_n(r-T)=\sum_{k}exp{-ikT}V_H^{k}(r)
!! Vc_n(r+T')=\sum_{k}exp{+ikT'}V_H^{k}(r)
!    do itloc=1,ntrloc
!      it=mpi_grid_map(ntr,dim_t,loc=itloc)
!      vtrc(:)=vtl(1,it)*avec(:,1)+vtl(2,it)*avec(:,2)+vtl(3,it)*avec(:,3)
!      expikt=exp(zi*dot_product(vkcnr(:,ik),vtrc(:)))
!      vcwanmt(:,:,:,itloc,n)=vcwanmt(:,:,:,itloc,n)+expikt*vckwanmt(:,:,:)
!      vcwanir(:,itloc,n)=vcwanir(:,itloc,n)+expikt*vckwanir(:)
!    enddo
!  enddo !n
!enddo !ikloc
!do n=1,nwann
!  do itloc=1,ntrloc 
!    call mpi_grid_reduce(vcwanmt(1,1,1,itloc,n),lmmaxvr*nrmtmax*natmtot,dims=(/dim_k/))
!    call mpi_grid_reduce(vcwanir(1,itloc,n),ngrtot,dims=(/dim_k/))
!  enddo
!enddo
!vcwanmt=vcwanmt/nkptnr
!vcwanir=vcwanir/nkptnr
!if (wproc) then
!  call timestamp(151,'done with Coulomb potential')
!endif


! for debug: sum over all WFs and translations to get the total Coulomb potential
!f2mt=zzero
!f2ir=zzero
!do itloc=1,ntrloc
!  do n=1,nwann
!    f2mt(:,:,:)=f2mt(:,:,:)+vcwanmt(:,:,:,itloc,n)
!    f2ir(:)=f2ir(:)+vcwanir(:,itloc,n)
!  enddo
!enddo
!do ir=1,nrmt(ias2is(1))
!  write(100,*)spr(ir,ias2is(1)),dreal(f2mt(1,ir,1)),dimag(f2mt(1,ir,1))
!enddo

allocate(rhowanmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nwann))
allocate(rhowanir(ngrtot,ntrloc,nwann))
! convert to spherical coordinates
do n=1,nwann
  do itloc=1,ntrloc
    call lfa_sht('B',wanmt(1,1,1,itloc,1,n),wanmt(1,1,1,itloc,1,n))
    rhowanmt(:,:,:,itloc,n)=dconjg(wanmt(:,:,:,itloc,1,n))*wanmt(:,:,:,itloc,1,n)
    rhowanir(:,itloc,n)=dconjg(wanir(:,itloc,1,n))*wanir(:,itloc,1,n)
  enddo
enddo

m=max(lmmaxvr,ngrtot)
allocate(vx(m),vc(m))
allocate(f3(m),f4(m))

! add XC potential to Coulomb
!do n=1,nwann
!  do itloc=1,ntrloc
!    do ias=1,natmtot
!      do ir=1,nrmt(ias2is(ias))
!        call xcifc(xctype,n=lmmaxvr,rho=dreal(rhowanmt(:,ir,ias,itloc,n)),ex=f3,ec=f4,vx=vx,vc=vc)
!        vcwanmt(:,ir,ias,itloc,n)=vcwanmt(:,ir,ias,itloc,n)+vc(1:lmmaxvr)+vx(1:lmmaxvr)
!      enddo
!    enddo
!    call xcifc(xctype,n=ngrtot,rho=dreal(rhowanir(:,itloc,n)),ex=f3,ec=f4,vx=vx,vc=vc)
!    vcwanir(:,itloc,n)=vcwanir(:,itloc,n)+vc(1:ngrtot)+vx(1:ngrtot)
!  enddo
!enddo
deallocate(vx,vc,f3,f4)
! compute bare Coulomb interaction
!allocate(ubare(nwann,nwann,ntr))
!ubare=zzero
!do it=1,ntr
!  do n1=1,nwann
!    do n2=1,nwann
!      ubare(n1,n2,it)=lfa_dotp(.false.,vtl(1,it),vcwanmt(1,1,1,1,n1),&
!        vcwanir(1,1,n1),rhowanmt(1,1,1,1,n2),rhowanir(1,1,n2))
!    enddo
!  enddo
!enddo
!if (wproc) then
!  write(151,*)
!  write(151,'("Bare Coulomb interaction")')
!  do it=1,ntr
!    write(151,'("  translation : ",3I4)')vtl(:,it)
!    write(151,'("    real part : ")')  
!    do n1=1,nwann
!      write(151,'(4X,255F12.6)')(dreal(ubare(n1,n2,it)),n2=1,nwann)
!    enddo
!    write(151,'("    image part : ")')  
!    do n1=1,nwann
!      write(151,'(4X,255F12.6)')(dimag(ubare(n1,n2,it)),n2=1,nwann)
!    enddo
!  enddo
!  write(151,*)
!endif
!if (wproc) then
!  call timestamp(151,'done with bare Coulomb integrals')
!endif


! multiply potential by Wannier function
do n=1,nwann
  do itloc=1,ntrloc
    vcwanmt(:,:,:,itloc,n)=vcwanmt(:,:,:,itloc,n)*wanmt(:,:,:,itloc,1,n)
    vcwanir(:,itloc,n)=vcwanir(:,itloc,n)*wanir(:,itloc,1,n)
  enddo
enddo

allocate(vsic(nwann,nwann,ntr))
vsic=zzero

! compute matrix elements of sic potential
! vsic = <w_n1|v_n1|w_{n2,T}>
do it=1,ntr
  do n1=1,nwann
    do n2=1,nwann
      vsic(n1,n2,it)=lfa_dotp(.false.,vtl(1,it),vcwanmt(1,1,1,1,n1),&
        vcwanir(1,1,n1),wanmt(1,1,1,1,1,n2),wanir(1,1,1,n2))
    enddo
  enddo
enddo
if (wproc) then
  write(151,*)
  write(151,'("SIC potential matrix elements <w_n1|v_n1|w_{n2,T}>")')
  do it=1,ntr
    write(151,'("  translation : ",3I4)')vtl(:,it)
    write(151,'("    real part : ")')  
    do n1=1,nwann
      write(151,'(4X,255F12.6)')(dreal(vsic(n1,n2,it)),n2=1,nwann)
    enddo
    write(151,'("    image part : ")')  
    do n1=1,nwann
      write(151,'(4X,255F12.6)')(dimag(vsic(n1,n2,it)),n2=1,nwann)
    enddo
  enddo
  write(151,*)
  write(151,'("SIC ""localization criterion"" <w_n1|v_n1-v_{n2,T}|w_{n2,T}>")')
  do it=1,ntr
    jt=ivtit(-vtl(1,it),-vtl(2,it),-vtl(3,it))
    write(151,'("  translation : ",3I4)')vtl(:,it)
    write(151,'("    real part : ")')  
    do n1=1,nwann
      write(151,'(4X,255F12.6)')(dreal(vsic(n1,n2,it)-dconjg(vsic(n2,n1,jt))),n2=1,nwann)
    enddo
    write(151,'("    image part : ")')  
    do n1=1,nwann
      write(151,'(4X,255F12.6)')(dimag(vsic(n1,n2,it)-dconjg(vsic(n2,n1,jt))),n2=1,nwann)
    enddo
  enddo
endif
if (wproc) then
  call timestamp(151,'done with SIC potential')
endif


!if (wproc) write(151,*)
!allocate(f1mt(lmmaxvr,nrmtmax,natmtot,ntrloc))
!allocate(f1ir(ngrtot,ntrloc))
!do n1=1,nwann
!  call lfa_prod(wanmt(1,1,1,1,1,n1),wanir(1,1,1,n1),wanmt(1,1,1,1,1,n1),&
!    wanir(1,1,1,n1),f1mt,f1ir)
!  z1=lfa_dotp(.false.,(/0,0,0/),identmt,identir,f1mt,f1ir)
!  if (wproc) write(151,'("norm : ",2G18.10)')z1
!enddo
!
!




!call charge
!
!
!allocate(f2mt(lmmaxvr,nrmtmax,natmtot))
!allocate(f2ir(ngrtot))
!f2mt=zzero
!f2ir=zzero
!!allocate(f1ir(ngrtot,ntrloc))
!do n1=1,nwann
!  call lfa_prod(wanmt(1,1,1,1,1,n1),wanir(1,1,1,n1),wanmt(1,1,1,1,1,n1),&
!    wanir(1,1,1,n1),f1mt,f1ir)
!  !z1=lfa_dotp(.false.,(/0,0,0/),identmt,identir,f1mt,f1ir)
!  !if (wproc) write(151,'("norm : ",2G18.10)')z1
!  do itrloc=1,ntrloc
!    f2mt(:,:,:)=f2mt(:,:,:)+f1mt(:,:,:,itrloc)
!    f2ir(:)=f2ir(:)+f1ir(:,itrloc)
!  enddo
!enddo
!call lfa_sht('F',f2mt,f2mt)
!rhomt=f2mt
!rhoir=f2ir
!call charge
!
!!if (wproc) then
!!  write(151,*)
!!  write(151,'("charge diff : ",G18.10)')sum(abs(rhomt-f2mt))
!!endif

if (wproc) close(151)
return
end

