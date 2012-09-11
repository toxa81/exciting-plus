subroutine sic_main
use modmain
use mod_nrkp
use mod_sic
use mod_wannier
use mod_hdf5
implicit none
integer n,j,ispn
integer ias,lm,iter
character*20 c1,c2,c3
character, parameter :: orbc(4)=(/'s','p','d','f'/)
character*2, parameter :: spinc(2)=(/'up','dn'/)
!
sic=.true.
! initialise universal variables
call init0
call init1
if (.not.mpi_grid_in()) return
wproc=mpi_grid_root()
! read the density and potentials from file
call readstate
! generate radial functions
call genradf
! generate wave-functions for all k-points in BZ
call genwfnr(-1,.false.)  
call elk_m_init
#ifdef _MAD_
!call madness_init_box
!call madness_gen_hpot(1)
#endif
call sic_create_hdf5
sic_write_rholm=.true.
sic_write_vlm=.true.
sic_write_vhlm=.true.
sic_write_vxclm=.true.
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
        &sic_wantran%iwan(j),wanpos(:,n)
    else
      ias=wan_info(wi_atom,n)
      lm=wan_info(wi_lm,n)
      ispn=wan_info(wi_spin,n)
      write(c1,'(I6)')ias2ia(ias)
      write(c2,'(I1)')lm2m(lm)+lm2l(lm)+1
      c3=trim(spsymb(ias2is(ias)))//trim(adjustl(c1))//"-"//&
        &orbc(lm2l(lm)+1)//trim(adjustl(c2))
      if (spinpol) c3=trim(adjustl(c3))//"-"//spinc(ispn)
      write(151,'("  j : ",I4,"    n : ",I4,"    ! ",A,"  at ",3G18.10)')j,&
        &sic_wantran%iwan(j),trim(c3),wanpos(:,n)
    endif
  enddo
  write(151,*)
  write(151,'("expansion radius for Wannier functions : ",F12.6)')s_rmax
  write(151,'("radius for radial integrals            : ",F12.6)')s_rmin
  write(151,'("number of radial points                : ",I6)')s_nr
  write(151,'("number of spherical points             : ",I6)')s_ntp
  write(151,'("maximum angular momentum               : ",I6)')lmaxwan
  write(151,'("cutoff radius for SIC matrix elements  : ",F12.6)')sic_me_cutoff
  write(151,'("number of Wannier transitions          : ",I6)')sic_wantran%nwt
  write(151,'("number of on-site Wannier transitions  : ",I6)')sic_wantran%nwt0
  !write(151,*)
  !write(151,'("LDA energies of Wannier functions")')
  !do j=1,sic_wantran%nwan
  !  n=sic_wantran%iwan(j)
  !  write(151,'("  n : ",I4,"    sic_wan_e0 : ",F12.6)')n,sic_wan_e0(n)
  !enddo
  call flushifc(151)
endif

do iter=1,sic_niter_umtrx
  if (iter.gt.1) then
    call sic_update_umtrx
  endif
  call wancnr_transform(sic_wan_umtrx)
  if (wproc) then
    write(151,*)
    write(151,'(80("="))')
    write(151,'("SIC localization step : ",I4)')iter
    write(151,'(80("="))')
    write(151,'("energies of Wannier functions")')
    do j=1,sic_wantran%nwan
      n=sic_wantran%iwan(j)
      write(151,'("  n : ",I4,"    wann_ene : ",F12.6)')n,wann_ene(n)
    enddo
  endif
! generate Wannier functions and corresponding potential
  call sic_wan_pot(151)
! matrix elements
  call sic_genvme(151,iter)
enddo
! write to HDF5 file
call sic_write_data
call mpi_grid_barrier
if (wproc) then
  call timestamp(151,"Done.")
  close(151)
endif
return
end

