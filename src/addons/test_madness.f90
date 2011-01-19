#ifdef _MAD_
subroutine test_madness
use modmain
use mod_nrkp
use mod_wannier
use mod_madness
implicit none
integer i,ik,ikloc,ig
real(8) d1,vrc(3),val(2)

call init0
call init1
write(*,*)"Hello from Elk from process ",iproc
! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! generate the core wavefunctions and densities
call gencore
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
call getufr
call genufrp
wproc=.false.
call genwfnr(-1,.false.)

if (allocated(m_ngknr)) deallocate(m_ngknr)
allocate(m_ngknr(nkptnr))
m_ngknr=0
if (allocated(m_igkignr)) deallocate(m_igkignr)
allocate(m_igkignr(ngkmax,nkptnr))
m_igkignr=0
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  m_ngknr(ik)=ngknr(ikloc)
  m_igkignr(:,ik)=igkignr(:,ikloc)
enddo
call mpi_grid_reduce(m_ngknr(1),nkptnr,dims=(/dim_k/),all=.true.)
call mpi_grid_reduce(m_igkignr(1,1),ngkmax*nkptnr,dims=(/dim_k/),all=.true.)
m_ngvec=0
do ik=1,nkptnr
  do ig=1,m_ngknr(ik)
    m_ngvec=max(m_ngvec,m_igkignr(ig,ik))
  enddo
enddo

if (allocated(m_wann_unkmt)) deallocate(m_wann_unkmt)
allocate(m_wann_unkmt(lmmaxvr,nufrmax,natmtot,nspinor,nkptnr))
if (allocated(m_wann_unkit)) deallocate(m_wann_unkit)
allocate(m_wann_unkit(ngkmax,nspinor,nkptnr))


!call timer_start(40,reset=.true.)
!vrc=0.d0
!call elk_load_wann_unk(1)
!do i=1,10000
!  vrc(1)=1.251d0
!  call elk_wan_rho(1,8.d0,vrc,val)
!enddo
!call timer_stop(40)
!write(*,*)timer_get_value(40)
!allocate(wann_unkmt1(lmmaxvr,nufrmax,natmtot,nspinor,nwantot))
!allocate(wann_unkit1(ngkmax,nspinor,nwantot))

call madness_init_box
!call mpi_world_barrier
call madness_resolve_wannier(nwantot)

!call madness_genpot
!call madness_getpot(d1)
!write(*,*)d1
call bstop
return
end
#endif
