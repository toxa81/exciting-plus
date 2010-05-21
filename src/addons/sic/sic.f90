subroutine sic
use modmain
use mod_lfa
use mod_nrkp
use modxcifc
implicit none

integer n
integer ik,ikloc,j,ig,sz,i
integer n1,ispn
integer itr,ntrloc,it,itloc,ir,m,ias
real(8) t1,t2
integer v1l(3)
character*12 c1,c2

real(8), allocatable :: vgq0c(:,:,:)
real(8), allocatable :: vhgq0(:,:)
real(8) a0

! arrays for Wannier functions
complex(8), allocatable :: wanmt(:,:,:,:,:,:)
complex(8), allocatable :: wanir(:,:,:,:)
complex(8), allocatable :: wanmt0(:,:,:,:,:)
complex(8), allocatable :: wanir0(:,:,:)
! Hartree potential of a Wanier state
complex(8), allocatable :: vhwanmt(:,:,:,:,:)
complex(8), allocatable :: vhwanir(:,:,:)
! plane-wave
complex(8), allocatable :: pwmt(:,:,:)
complex(8), allocatable :: pwir(:)

complex(8), allocatable :: megqwan1(:,:,:) 

complex(8), allocatable :: rhowanmt(:,:,:)
complex(8), allocatable :: rhowanir(:)
complex(8), allocatable :: vsic(:)
real(8), allocatable :: f3(:),f4(:)

integer nvq0loc,iqloc,iq
real(8), allocatable :: vx(:),vc(:)
complex(8) z1
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

call lfa_init(1)
call genwfnr(151,.true.)  
call init_qbz(.true.)
call getngvecme((/0,0,0/))

wproc=mpi_grid_root()
if (wproc) then
  open(151,file='SIC.OUT',form='FORMATTED',status='REPLACE')
endif
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
if (wproc) then
  call timestamp(151,'done with wavefunctions')
  call flushifc(151)
endif
! get all Wannier transitions
call getimegqwan(.true.)
! create q-directories
!if (mpi_grid_root()) then
!  call system("mkdir -p q")
!  do iq=1,nvq0
!    call getqdir(iq,ivq0m_list(:,iq),qnm)
!    call system("mkdir -p "//trim(qnm))
!  enddo
!endif
wannier_megq=.true.
all_wan_ibt=.true.
! distribute q-vectors along 3-rd dimention
nvq0loc=mpi_grid_map(nvq0,dim_q)
! distribute translations along 3-rd dimention
ntrloc=mpi_grid_map(ntr,dim_t)
allocate(megqwan1(nwann,ngvecme,nvq0))
! loop over q-points
do iqloc=1,nvq0loc
  iq=mpi_grid_map(nvq0,dim_q,loc=iqloc)
  call genmegq(iq,.false.)
! save <n,T=0|e^{-i(G+q)r|n,T=0>
  do n=1,nwann
    megqwan1(n,:,iq)=megqwan(idxmegqwan(n,n,0,0,0),:)
  enddo
enddo
if (wproc) then
  call timestamp(151,'done with q-vectors')
endif
call mpi_grid_reduce(megqwan1(1,1,1),nwann*ngvecme*nvq0,dims=(/dim_q/), &
  all=.true.)
! deallocate unnecessary wave-functions
deallocate(wfsvmtloc)
deallocate(wfsvitloc)
deallocate(evecfvloc)
deallocate(evecsvloc)
deallocate(wann_c)
! restore wproc
wproc=mpi_grid_root()
! generate G+q vectors   
allocate(vgq0c(3,ngvecme,nvq0))
allocate(vhgq0(ngvecme,nvq0))
do iq=1,nvq0
  vq0l(:)=1.d0*(ivq0m_list(:,iq))/ngridk(:)+1d-12
  call r3mv(bvec,vq0l,vq0c)
  do ig=1,ngvecme
    if (ivq0m_list(1,iq).eq.0.and.&
        ivq0m_list(2,iq).eq.0.and.&
        ivq0m_list(3,iq).eq.0) then
      if (ig.eq.1) then
        vgq0c(:,ig,iq)=q0gamma(:,iq)
        a0=a0gamma(iq)
      else
        vgq0c(:,ig,iq)=vgc(:,ig)
        a0=0.125d0
      endif
    else
      vgq0c(:,ig,iq)=vgc(:,ig)+vq0c(:)
      a0=1.d0
    endif
! setup 4Pi/|G+q|^2 array
    vhgq0(ig,iq)=a0*fourpi/dot_product(vgq0c(:,ig,iq),vgq0c(:,ig,iq))
  enddo !ig
enddo !iq
! allocate arrays for plane-wave
allocate(pwmt(lmmaxvr,nrmtmax,natmtot))
allocate(pwir(ngrtot))
! generate Hartree potential
allocate(vhwanmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nwann))
allocate(vhwanir(ngrtot,ntrloc,nwann))
vhwanmt=zzero
vhwanir=zzero
do itloc=1,ntrloc
  itr=mpi_grid_map(ntr,dim_t,loc=itloc)
  do ig=1,ngvecme
    do iq=1,nvq0
      call genpw(vtl(1,itr),vgq0c(1,ig,iq),pwmt,pwir)
      do n=1,nwann
        vhwanmt(:,:,:,itloc,n)=vhwanmt(:,:,:,itloc,n)+&
          megqwan1(n,ig,iq)*vhgq0(ig,iq)*pwmt(:,:,:)
        vhwanir(:,itloc,n)=vhwanir(:,itloc,n)+&
          megqwan1(n,ig,iq)*vhgq0(ig,iq)*pwir(:)
      enddo
    enddo
  enddo
enddo
vhwanmt=vhwanmt/nkptnr/omega
vhwanir=vhwanir/nkptnr/omega
if (wproc) then
  call timestamp(151,'done with Hartree potential')
endif
deallocate(vgq0c,vhgq0,pwmt,pwir,megqwan1)
! generate Wannier functions on a mesh
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
! check orthonormality
t1=0.d0
t2=0.d0
do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  v1l(:)=imegqwan(3:5,i)
  z1=0
  do ispn=1,nspinor
    z1=z1+lfa_dotp(.true.,v1l,wanmt(1,1,1,1,ispn,n),wanir(1,1,ispn,n),&
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


! test
!allocate(pwmt1(lmmaxvr,nrmtmax,natmtot,ntrloc))
!allocate(pwir1(ngrtot,ntrloc))
!do n=1,nwann
!  do itloc=1,ntrloc
!    call lfa_sht('B',wanmt(1,1,1,itloc,1,n),wanmt(1,1,1,itloc,1,n))
!  enddo
!enddo
!do ig=1,ngvecme
!  do iq=1,nvq0
!    do itloc=1,ntrloc
!      itr=mpi_grid_map(ntr,dim_t,loc=itloc)
!      call genpw(vtl(1,itr),vgq0c(1,ig,iq),pwmt1(1,1,1,itloc),pwir1(1,itloc))
!      call lfa_sht('B',pwmt1(1,1,1,itloc),pwmt1(1,1,1,itloc))
!      pwmt1(:,:,:,itloc)=pwmt1(:,:,:,itloc)*wanmt(:,:,:,itloc,1,1)
!      pwir1(:,itloc)=pwir1(:,itloc)*wanir(:,itloc,1,1)    
!    enddo
!    zt1=lfa_dotp(.false.,(/0,0,0/),pwmt1,pwir1,wanmt(1,1,1,1,1,1),wanir(1,1,1,1))
!    if (wproc) write(*,*)'ig=',ig,'iq=',iq,'  ',zt1
!  enddo
!enddo
!call mpi_grid_barrier
!call pstop


! convert to spherical coordinates
do n=1,nwann
  do itloc=1,ntrloc
    call lfa_sht('B',wanmt(1,1,1,itloc,1,n),wanmt(1,1,1,itloc,1,n))
    call lfa_sht('B',vhwanmt(1,1,1,itloc,n),vhwanmt(1,1,1,itloc,n))  
  enddo
enddo

m=max(lmmaxvr,ngrtot)
allocate(rhowanmt(lmmaxvr,nrmtmax,natmtot))
allocate(rhowanir(ngrtot))
allocate(vx(m),vc(m))
allocate(f3(m),f4(m))
! add XC potential to Coulomb
do n=1,nwann
  do itloc=1,ntrloc
    rhowanmt(:,:,:)=dconjg(wanmt(:,:,:,itloc,1,n))*wanmt(:,:,:,itloc,1,n)
    rhowanir(:)=dconjg(wanir(:,itloc,1,n))*wanir(:,itloc,1,n)
    do ias=1,natmtot
      do ir=1,nrmt(ias2is(ias))
        call xcifc(xctype,n=lmmaxvr,rho=dreal(rhowanmt(:,ir,ias)),&
          ex=f3,ec=f4,vx=vx,vc=vc)
        vhwanmt(:,ir,ias,itloc,n)=vhwanmt(:,ir,ias,itloc,n)+vc(1:lmmaxvr)+&
          vx(1:lmmaxvr)
      enddo
    enddo
    call xcifc(xctype,n=ngrtot,rho=dreal(rhowanir(:)),ex=f3,ec=f4,&
      vx=vx,vc=vc)
    vhwanir(:,itloc,n)=vhwanir(:,itloc,n)+vc(1:ngrtot)+vx(1:ngrtot)
  enddo
enddo
deallocate(vx,vc,f3,f4,rhowanmt,rhowanir)

! multiply potential by Wannier function and change sign
do n=1,nwann
  do itloc=1,ntrloc
    vhwanmt(:,:,:,itloc,n)=-vhwanmt(:,:,:,itloc,n)*wanmt(:,:,:,itloc,1,n)
    vhwanir(:,itloc,n)=-vhwanir(:,itloc,n)*wanir(:,itloc,1,n)
  enddo
enddo

allocate(vsic(nmegqwan))
! compute matrix elements of SIC potential
! vsic = <w_n|v_n|w_{n1,T}>
do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  v1l(:)=imegqwan(3:5,i)
  vsic(i)=lfa_dotp(.false.,v1l,vhwanmt(1,1,1,1,n),vhwanir(1,1,n),&
    wanmt(1,1,1,1,1,n1),wanir(1,1,1,n1))
enddo
if (wproc) then
  write(151,*)
  write(151,'("Number of Wannier transitions : ",I6)')nmegqwan
  write(151,'("Matrix elements of SIC potential \&
    &(n n1  T  <w_n|v_n|w_{n1,T}>)")')
  do i=1,nmegqwan
    write(151,'(I4,4X,I4,4X,3I3,4X,2G18.10)')imegqwan(:,i),&
      dreal(vsic(i)),dimag(vsic(i))
  enddo
endif

t2=0.d0
do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  v1l(:)=imegqwan(3:5,i)
  t2=max(t2,abs(vsic(i)-dconjg(vsic(idxmegqwan(n1,n,-v1l(1),-v1l(2),-v1l(3))))))
enddo
if (wproc) then
  write(151,*)
  write(151,'("Maximum deviation from ""localization criterion"" : ",F12.6)')t2
  call flushifc(151)
endif
if (wproc) then
  call hdf5_create_file("sic.hdf5")
  call hdf5_create_group("sic.hdf5","/","vwan")
  do n=1,nwann
    write(c1,'("n",I4.4)')n
    call hdf5_create_group("sic.hdf5","/vwan",trim(adjustl(c1)))
    do it=1,ntr
      write(c2,'("t",I4.4)')it
      call hdf5_create_group("sic.hdf5","/vwan/"//trim(adjustl(c1)),&
        trim(adjustl(c2)))
    enddo
  enddo
  call hdf5_write("sic.hdf5","/","nmegqwan",nmegqwan)
  call hdf5_write("sic.hdf5","/","imegqwan",imegqwan(1,1),(/5,nmegqwan/))
  call hdf5_write("sic.hdf5","/","vsic",vsic(1),(/nmegqwan/))
endif
if (mpi_grid_side(dims=(/dim_t/))) then
  do i=0,mpi_grid_size(dim_t)-1
    if (mpi_grid_x(dim_t).eq.i) then
      do itloc=1,ntrloc
        itr=mpi_grid_map(ntr,dim_t,loc=itloc)
        do n=1,nwann
          write(c1,'("n",I4.4)')n
          write(c2,'("t",I4.4)')itr
          call hdf5_write("sic.hdf5","/vwan/"//trim(adjustl(c1))//"/"//&
            trim(adjustl(c2)),"vwanmt",vhwanmt(1,1,1,itloc,n),&
            (/lmmaxvr,nrmtmax,natmtot/))
          call hdf5_write("sic.hdf5","/vwan/"//trim(adjustl(c1))//"/"//&
            trim(adjustl(c2)),"vwanir",vhwanir(1,itloc,n),&
            (/ngrtot/))
        enddo
      enddo
    endif
    call mpi_grid_barrier(dims=(/dim_t/))
  enddo
endif
          
          
if (wproc) close(151)
return
end

