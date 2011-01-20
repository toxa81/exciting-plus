#ifdef _MAD_
subroutine test_madness
use modmain
use mod_nrkp
use mod_wannier
use mod_madness
implicit none
integer i,ik,ikloc,ig
real(8) d1,d2,vrc(3)

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

!
!
!vrc=0.d0
!call elk_load_wann_unk(1)
!do i=1,10000
!  vrc(1)=15d0*i/10000.d0
!  call elk_wan_rho(1,10.d0,vrc,d1)
!  write(99,*)vrc(1),d1
!enddo


call madness_init_box

call elk_load_wann_unk(1)
call madness_gen_hpot(1)
open(90,file="wann_s.dat",form="formatted",status="replace")
vrc=0.d0
do i=1,10000
  vrc(1)=10d0*i/10000.d0
  call madness_get_hpot(vrc,d1)
  call madness_get_rho(vrc,d2)
  write(90,'(3G18.10)')vrc(1),d1,d2
enddo
close(90)

call elk_load_wann_unk(4)
call madness_gen_hpot(4)
open(90,file="wann_px.dat",form="formatted",status="replace")
vrc=0.d0
do i=1,10000
  vrc(1)=10d0*i/10000.d0
  call madness_get_hpot(vrc,d1)
  call madness_get_rho(vrc,d2)
  write(90,'(3G18.10)')vrc(1),d1,d2
enddo
close(90)

call bstop
return
end
#endif
