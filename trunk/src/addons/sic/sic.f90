subroutine sic
use modmain
implicit none

integer i1,i2,i3,n
integer ik,ikloc,ik1,j,ig,sz,i,isym,jtrloc
integer ispn,n1,n2
complex(8) znorm
integer ntr,itr,ntrloc,itrloc
integer, allocatable :: vtrl(:,:)
integer, parameter :: trmax=0

complex(8), allocatable :: wanmt(:,:,:,:,:,:)
complex(8), allocatable :: wanir(:,:,:,:)
complex(8), allocatable :: wanmt0(:,:,:,:,:)
complex(8), allocatable :: wanir0(:,:,:)

integer, allocatable :: ngknr(:)
integer, allocatable :: igkignr(:,:)
real(8), allocatable :: vgklnr(:,:,:)
real(8), allocatable :: vgkcnr(:,:,:)
real(8), allocatable :: gknr(:,:)
real(8), allocatable :: tpgknr(:,:,:)
complex(8), allocatable :: sfacgknr(:,:,:)
complex(8), allocatable :: ylmgknr(:,:,:)

complex(8), allocatable :: wfsvmtloc(:,:,:,:,:,:)
complex(8), allocatable :: wfsvitloc(:,:,:,:)
complex(8), allocatable :: wfsvcgloc(:,:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: tp(:,:)
complex(8), allocatable :: ylmtp(:)
real(8) vtrc(3)
complex(8), allocatable :: ovlm(:,:)
complex(8), external :: zfinp_





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

! get energies of states in reduced part of BZ
call timer_start(3,reset=.true.)
if (wproc) then
  write(151,*)
  write(151,'("Reading energies of states")')
  call flushifc(151)
! read from IBZ
  do ik=1,nkpt
    call getevalsv(vkl(1,ik),evalsv(1,ik))
  enddo
endif
call mpi_grid_bcast(evalsv(1,1),nstsv*nkpt)
allocate(lr_evalsvnr(nstsv,nkptnr))
lr_evalsvnr=0.d0
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  call findkpt(vklnr(1,ik),isym,ik1) 
  lr_evalsvnr(:,ik)=evalsv(:,ik1)
enddo
call timer_stop(3)
if (wproc) then
  write(151,'("Done in ",F8.2," seconds")')timer_get_value(3)
  call timestamp(151)
  call flushifc(151)
endif

! generate G+k vectors for entire BZ (this is required to compute 
!   wave-functions at each k-point)
allocate(vgklnr(3,ngkmax,nkptnrloc))
allocate(vgkcnr(3,ngkmax,nkptnrloc))
allocate(gknr(ngkmax,nkptnrloc))
allocate(tpgknr(2,ngkmax,nkptnrloc))
allocate(ngknr(nkptnrloc))
allocate(sfacgknr(ngkmax,natmtot,nkptnrloc))
allocate(igkignr(ngkmax,nkptnrloc))
allocate(ylmgknr(lmmaxvr,ngkmax,nkptnrloc))
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  call gengpvec(vklnr(1,ik),vkcnr(1,ik),ngknr(ikloc),igkignr(1,ikloc), &
    vgklnr(1,1,ikloc),vgkcnr(1,1,ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc))
  call gensfacgp(ngknr(ikloc),vgkcnr(1,1,ikloc),ngkmax,sfacgknr(1,1,ikloc))
  do ig=1,ngknr(ikloc)
    call genylm(lmaxvr,tpgknr(1,ig,ikloc),ylmgknr(1,ig,ikloc))
  enddo
enddo
allocate(wfsvmtloc(lmmaxvr,nrfmax,natmtot,nspinor,nstsv,nkptnrloc))
allocate(wfsvitloc(ngkmax,nspinor,nstsv,nkptnrloc))
allocate(wfsvcgloc(ngkmax,nspinor,nstsv,nkptnrloc))  
allocate(evecfvloc(nmatmax,nstfv,nspnfv,nkptnrloc))
allocate(evecsvloc(nstsv,nstsv,nkptnrloc))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
if (wproc) then
  sz=lmmaxvr*nrfmax*natmtot*nstsv*nspinor
  sz=sz+ngkmax*nstsv*nspinor
  sz=sz+nmatmax*nstfv*nspnfv
  sz=sz+nstsv*nstsv
  sz=16*sz*nkptnrloc/1024/1024
  write(151,*)
  write(151,'("Size of wave-function arrays (MB) : ",I6)')sz
  write(151,*)
  write(151,'("Reading eigen-vectors")')
  call flushifc(151)
endif
call timer_start(1,reset=.true.)
! read eigen-vectors
if (mpi_grid_side(dims=(/dim_k/))) then
  do i=0,mpi_grid_size(dim_k)-1
    if (i.eq.mpi_grid_x(dim_k)) then
      do ikloc=1,nkptnrloc
        ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
        call getevecfv(vklnr(1,ik),vgklnr(1,1,ikloc),evecfvloc(1,1,1,ikloc))
        call getevecsv(vklnr(1,ik),evecsvloc(1,1,ikloc))
      enddo !ikloc
    endif
    if (.not.parallel_read) call mpi_grid_barrier(dims=(/dim_k/))
  enddo
endif !mpi_grid_side(dims=(/dim_k/)
call mpi_grid_barrier
!call mpi_grid_bcast(evecfvloc(1,1,1,1),nmatmax*nstfv*nspnfv*nkptnrloc,&
!  dims=(/dim2,dim3/))
!call mpi_grid_bcast(evecsvloc(1,1,1),nstsv*nstsv*nkptnrloc,&
!  dims=(/dim2,dim3/))
! transform eigen-vectors
wfsvmtloc=zzero
wfsvitloc=zzero
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
! get apw coeffs 
  call match(ngknr(ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc),        &
    sfacgknr(1,1,ikloc),apwalm)
! generate wave functions in muffin-tins
  call genwfsvmt(lmaxvr,lmmaxvr,ngknr(ikloc),evecfvloc(1,1,1,ikloc), &
    evecsvloc(1,1,ikloc),apwalm,wfsvmtloc(1,1,1,1,1,ikloc))
! generate wave functions in interstitial
  call genwfsvit(ngknr(ikloc),evecfvloc(1,1,1,ikloc), &
    evecsvloc(1,1,ikloc),wfsvitloc(1,1,1,ikloc))
enddo !ikloc
call timer_stop(1)
if (wproc) then
  write(151,'("Done in ",F8.2," seconds")')timer_get_value(1)
  call timestamp(151)
  call flushifc(151)
endif
! generate Wannier function expansion coefficients
if (wannier) then
  call timer_start(1,reset=.true.)
  if (allocated(wann_c)) deallocate(wann_c)
  allocate(wann_c(nwann,nstsv,nkptnrloc))
  if (allocated(wann_unkmt)) deallocate(wann_unkmt)
  allocate(wann_unkmt(lmmaxvr,nrfmax,natmtot,nspinor,nwann,nkptnrloc))
  if (allocated(wann_unkit)) deallocate(wann_unkit)
  allocate(wann_unkit(ngkmax,nspinor,nwann,nkptnrloc))
  wann_unkmt=zzero
  wann_unkit=zzero
  if (wproc) then
    write(151,*)
    write(151,'("Generating Wannier functions")')
    call flushifc(151)
  endif !wproc1
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    call genwann_c(ik,lr_evalsvnr(1,ik),wfsvmtloc(1,1,1,1,1,ikloc),&
      wann_c(1,1,ikloc))  
    do n=1,nwann
      do j=1,nstsv
        wann_unkmt(:,:,:,:,n,ikloc)=wann_unkmt(:,:,:,:,n,ikloc) + &
          wfsvmtloc(:,:,:,:,j,ikloc)*wann_c(n,j,ikloc)
        wann_unkit(:,:,n,ikloc)=wann_unkit(:,:,n,ikloc) + &
          wfsvitloc(:,:,j,ikloc)*wann_c(n,j,ikloc)
      enddo
    enddo
  enddo !ikloc
endif !wannier
! after optinal band disentanglement we can finally synchronize all eigen-values
!   and compute band occupation numbers 
call mpi_grid_reduce(lr_evalsvnr(1,1),nstsv*nkptnr,dims=(/dim_k/),all=.true.)
allocate(lr_occsvnr(nstsv,nkptnr))
call occupy2(nkptnr,wkptnr,lr_evalsvnr,lr_occsvnr)

deallocate(wfsvmtloc)
deallocate(wfsvitloc)
deallocate(wfsvcgloc)  
deallocate(evecfvloc)
deallocate(evecsvloc)
deallocate(apwalm)
deallocate(wann_c)


ntr=(2*trmax+1)**3
allocate(vtrl(3,ntr))
n=0
do i1=-trmax,trmax
  do i2=-trmax,trmax
    do i3=-trmax,trmax
      n=n+1
      vtrl(:,n)=(/i1,i2,i3/)
    enddo
  enddo
enddo
ntrloc=mpi_grid_map(ntr,dim_k)
if (wproc) then
  write(151,*)
  write(151,'("Number of translations : ",I4)')ntr
endif

allocate(wanmt(lmmaxvr,nrmtmax,natmtot,nspinor,nwann,ntrloc))
allocate(wanir(ngrtot,nspinor,nwann,ntrloc))
allocate(wanmt0(lmmaxvr,nrmtmax,natmtot,nspinor,nwann))
allocate(wanir0(ngrtot,nspinor,nwann))
allocate(tp(2,lmmaxvr))
allocate(ylmtp(lmmaxvr))
call sphcover(lmmaxvr,tp)

do itr=1,ntr
  jtrloc=mpi_grid_map(ntr,dim_k,glob=itr,x=j)
  vtrc(:)=vtrl(1,itr)*avec(:,1)+vtrl(2,itr)*avec(:,2)+vtrl(3,itr)*avec(:,3)
  call gen_wann_func(vtrc,ngknr,vgkcnr,wanmt0,wanir0)
  call mpi_grid_reduce(wanmt0(1,1,1,1,1),lmmaxvr*nrmtmax*natmtot*nspinor*nwann,&
    dims=(/dim_k/),root=(/j/))
  call mpi_grid_reduce(wanir0(1,1,1),ngrtot*nspinor*nwann,&
    dims=(/dim_k/),root=(/j/))
  if (mpi_grid_x(dim_k).eq.j) then
    wanmt(:,:,:,:,:,jtrloc)=wanmt0(:,:,:,:,:)
    wanir(:,:,:,jtrloc)=wanir0(:,:,:)
  endif
enddo !itrloc

allocate(ovlm(nwann,nwann))
ovlm=zzero
! calculate overlap matrix
do n1=1,nwann
  do n2=1,nwann
    znorm=zzero
    do itrloc=1,ntrloc
      do ispn=1,nspinor
        ovlm(n1,n2)=ovlm(n1,n2)+zfinp_(.true.,wanmt(1,1,1,ispn,n1,itrloc),&
          wanmt(1,1,1,ispn,n2,itrloc),wanir(1,ispn,n1,itrloc),&
          wanir(1,ispn,n2,itrloc))
      enddo
    enddo
  enddo
enddo
call mpi_grid_reduce(ovlm(1,1),nwann*nwann,dims=(/dim_k/))
if (wproc) then
  write(151,*)
  write(151,'("Overlap matrix")')
  write(151,'("  Real part : ")')  
  do n1=1,nwann
    write(151,'(2X,255F12.6)')(dreal(ovlm(n1,n2)),n2=1,nwann)
  enddo
  write(151,'("  Image part : ")')  
  do n1=1,nwann
    write(151,'(2X,255F12.6)')(dimag(ovlm(n1,n2)),n2=1,nwann)
  enddo
endif

if (wproc) close(151)
return
end



subroutine gen_wann_func(vtrc,ngknr,vgkcnr,wanmt,wanir)
use modmain
implicit none
real(8), intent(in) :: vtrc(3)
integer, intent(in) :: ngknr(nkptnrloc)
real(8), intent(in) :: vgkcnr(3,ngkmax,nkptnrloc)
complex(8), intent(out) :: wanmt(lmmaxvr,nrmtmax,natmtot,nspinor,nwann)
complex(8), intent(out) :: wanir(ngrtot,nspinor,nwann)
integer ia,is,ias,ir,ir0,i1,i2,i3,ig,ikloc
integer io,lm,n,ispn,itmp(3)
real(8) v2(3),v3(3),r0,vr0(3)
complex(8) expikr
logical, external :: vrinmt
wanmt=zzero
wanir=zzero
! muffin-tin part
call timer_start(1,reset=.true.)
do ias=1,natmtot
  is=ias2is(ias)
  do ikloc=1,nkptnrloc
    expikr=exp(zi*dot_product(vkcnr(:,ikloc),vtrc(:)))/nkptnr
    do ir=1,nrmt(is)
      do lm=1,lmmaxvr
        do io=1,nrfl(lm2l(lm),is)
          do ispn=1,nspinor
            do n=1,nwann
              wanmt(lm,ir,ias,ispn,n)=wanmt(lm,ir,ias,ispn,n)+&
                expikr*wann_unkmt(lm,io,ias,ispn,n,ikloc)*&
                urf(ir,lm2l(lm),io,ias)
            enddo !n
          enddo !ispn
        enddo !io
      enddo !lm
    enddo !ir
  enddo !ikloc
enddo !ias 
call timer_stop(1)
! interstitial part
call timer_start(2,reset=.true.)
ir=0
do i3=0,ngrid(3)-1
  v2(3)=dble(i3)/dble(ngrid(3))
  do i2=0,ngrid(2)-1
    v2(2)=dble(i2)/dble(ngrid(2))
    do i1=0,ngrid(1)-1
      v2(1)=dble(i1)/dble(ngrid(1))
      ir=ir+1
      call r3mv(avec,v2,v3)
      if (.not.vrinmt(v3,is,ia,itmp,vr0,ir0,r0)) then
        v3(:)=v3(:)+vtrc(:)
        do ikloc=1,nkptnrloc
          do ig=1,ngknr(ikloc)
            expikr=exp(zi*dot_product(vgkcnr(:,ig,ikloc),v3(:)))
            do ispn=1,nspinor
              do n=1,nwann
                wanir(ir,ispn,n)=wanir(ir,ispn,n)+&
                  expikr*wann_unkit(ig,ispn,n,ikloc)
              enddo !in
            enddo !ispn
          enddo !ig
        enddo !ikloc
      endif
    enddo !i1
  enddo !i2
enddo !i3
wanir(:,:,:)=wanir(:,:,:)/sqrt(omega)/nkptnr
call timer_stop(2)
return
end





complex(8) function zfinp_(tsh,zfmt1,zfmt2,zfir1,zfir2)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tsh   : .true. if the muffin-tin functions are in spherical harmonics
!           (in,logical)
!   zfmt1 : first complex function in spherical harmonics/coordinates for all
!           muffin-tins (in,complex(lmmaxvr,nrcmtmax,natmtot))
!   zfmt2 : second complex function in spherical harmonics/coordinates for all
!           muffin-tins (in,complex(lmmaxvr,nrcmtmax,natmtot))
!   zfir1 : first complex interstitial function in real-space
!           (in,complex(ngrtot))
!   zfir2 : second complex interstitial function in real-space
!           (in,complex(ngrtot))
! !DESCRIPTION:
!   Calculates the inner product of two complex fuctions over the entire unit
!   cell. The muffin-tin functions should be stored on the coarse radial grid
!   and have angular momentum cut-off {\tt lmaxvr}. In the intersitial region,
!   the integrand is multiplied with the characteristic function, to remove the
!   contribution from the muffin-tin. See routines {\tt zfmtinp} and
!   {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created July 2004 (Sharma)
!EOP
!BOC
implicit none
! arguments
logical, intent(in) :: tsh
complex(8), intent(in) :: zfmt1(lmmaxvr,nrmtmax,natmtot)
complex(8), intent(in) :: zfmt2(lmmaxvr,nrmtmax,natmtot)
complex(8), intent(in) :: zfir1(ngrtot)
complex(8), intent(in) :: zfir2(ngrtot)
! local variables
integer is,ia,ias,ir
complex(8) zsum
! external functions
complex(8) zfmtinp
external zfmtinp
zsum=0.d0
! interstitial contribution
do ir=1,ngrtot
  zsum=zsum+cfunir(ir)*conjg(zfir1(ir))*zfir2(ir)
end do
zsum=zsum*omega/dble(ngrtot)
! muffin-tin contribution
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    zsum=zsum+zfmtinp(tsh,lmaxvr,nrmt(is),spr(:,is),lmmaxvr,zfmt1(:,:,ias), &
     zfmt2(:,:,ias))
  end do
end do
zfinp_=zsum
return
end function
