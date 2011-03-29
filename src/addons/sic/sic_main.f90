subroutine sic_main
use modmain
use mod_nrkp
use mod_sic
use mod_wannier
use mod_madness
implicit none
integer n,j,ik,ispn,ikloc,ig
integer ias,lm
real(8), allocatable :: laplsv(:) 
character*20 c1,c2,c3
character, parameter :: orbc(4)=(/'s','p','d','f'/)
character*2, parameter :: spinc(2)=(/'up','dn'/)
!
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
  write(151,*)
  write(151,'("number of included Wannier functions : ",I4)')sic_wantran%nwan
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    if (wannier_lc) then
      write(151,'("  j : ",I4,"    n : ",I4,"    !  at ",3G18.10)')j,&
        sic_wantran%iwan(j),wanpos(:,n)
    else
      ias=wan_info(1,n)
      lm=wan_info(2,n)
      ispn=wan_info(3,n)
      write(c1,'(I6)')ias2ia(ias)
      write(c2,'(I1)')lm2m(lm)+lm2l(lm)+1
      c3=trim(spsymb(ias2is(ias)))//trim(adjustl(c1))//"-"//&
        orbc(lm2l(lm)+1)//trim(adjustl(c2))
      if (spinpol) c3=trim(adjustl(c3))//"-"//spinc(ispn)
      write(151,'("  j : ",I4,"    n : ",I4,"    ! ",A,"  at ",3G18.10)')j,&
        sic_wantran%iwan(j),trim(c3),wanpos(:,n)
    endif
  enddo
  write(151,*)
  write(151,'("cutoff radius for Wannier functions   : ",F12.6)')sic_wan_cutoff
  write(151,'("number of radial points               : ",I6)')s_nr
  write(151,'("number of spherical points            : ",I6)')s_ntp
  write(151,'("maximum angular momentum              : ",I6)')lmaxwan
  write(151,'("cutoff radius for SIC matrix elements : ",F12.6)')sic_me_cutoff
  write(151,'("number of Wannier transitions         : ",I6)')sic_wantran%nwt
  write(151,'("number of translations for Bloch sums : ",I6)')sic_orbitals%ntr
  write(151,*)
  write(151,'("LDA energies of Wannier functions")')
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    write(151,'("  n : ",I4,"    sic_wann_e0 : ",F12.6)')n,sic_wann_e0(n)
  enddo
  call flushifc(151)
endif
if (wproc) then
  write(151,*)
  write(151,'(80("="))')
  write(151,'("generating wave-functions for all k-points")')
  write(151,'(80("="))')
endif
! generate wave-functions for all k-points in BZ
call genwfnr(151,.false.)  
! compute kinetic energy
!allocate(laplsv(nstsv))
!ekin_=0.d0
!do ikloc=1,nkptnrloc
!  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
!  call genlapl(ngknr(ikloc),igkignr(1,ikloc),vgkcnr(1,1,ikloc),&
!    wfsvmtnrloc(1,1,1,1,1,ikloc),wfsvitnrloc(1,1,1,ikloc),laplsv)
!  do j=1,nstsv
!    ekin_=ekin_+0.5d0*laplsv(j)*occsvnr(j,ik)
!  enddo
!enddo
!call mpi_grid_reduce(ekin_,dims=(/dim_k/))
!ekin_=ekin_/nkptnr
!if (wproc) then
!  write(151,*)
!  etot_=ekin_+engykncr+0.5d0*engyvcl+engyx+engyc+engymad-sic_energy_pot
!  write(151,'("kinetic energy : ",G18.10)')ekin_+engykncr
!  write(151,'("potential      : ",G18.10)')0.5d0*engyvcl+engyx+engyc
!  write(151,'("madelung       : ",G18.10)')engymad
!  write(151,'("sic correction : ",G18.10)')sic_energy_pot
!  write(151,'("total energy   : ",G18.10)')etot_
!  open(152,file="SIC_ETOT.OUT",form="FORMATTED",status="REPLACE")
!  write(152,'(G18.10)')etot_
!  close(152)
!endif

! get maximum number of G-vectors (this is required for the faster generation of 
!  Wannier functions in the interstitial region)
s_ngvec=0
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  do ig=1,ngknr(ikloc)
    s_ngvec=max(s_ngvec,igkignr(ig,ikloc))
  enddo
enddo

!call sic_readvwan
!call sic_test_fvprj(151)
!call sic_test_blochsum(1,.true.,"sic_blochsum_wan.out")

! init Madness related variables 
#ifdef _MAD_
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
call madness_init_box
#endif
! generate Wannier functions and corresponding potential
call sic_wan(151)
! matrix elements
call sic_me(151)
!call sic_test_blochsum(2,.false.,"sic_blochsum_wan_backward.out")
!call sic_diff_blochsum_mt
! write to HDF5 file
call sic_writevwan
if (wproc) then
  call timestamp(151,"Done.")
  close(151)
endif
return
end

