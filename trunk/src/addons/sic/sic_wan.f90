subroutine sic_wan(fout)
use modmain
use mod_sic
use mod_nrkp
use mod_wannier
use mod_linresp
implicit none
! arguments
integer, intent(in) :: fout
! local variables
integer n,ispn,vl(3),n1,i,j,j1,itp,ir,nwtloc,iloc,nrloc,irloc
real(8) t1,vrc(3),x(3),x2,pos1(3),pos2(3)
real(8) sic_epot_h,sic_epot_xc
complex(8) z1,wanval(nspinor)
real(8), allocatable :: wanprop(:,:)
complex(8), allocatable :: ovlp(:)
complex(8), allocatable :: wantp(:,:,:)
integer, parameter :: iovlp=1 
integer, parameter :: irms=1
complex(8), external :: zdotc
!
allocate(wantp(s_ntp,s_nr,nspinor))
allocate(wanprop(nwanprop,sic_wantran%nwan))
if (wproc) then
  write(fout,*)
  write(fout,'(80("="))')
  write(fout,'("generating Wannier functions on a spherical grid")')
  write(fout,'(80("="))')
endif
call timer_start(t_sic_wan,reset=.true.)
call timer_reset(t_sic_wan_gen)
call timer_reset(t_sic_wan_rms)
call timer_reset(t_sic_wan_pot)
s_wanlm=zzero
s_wvlm=zzero
wanprop=0.d0
nrloc=mpi_grid_map(s_nr,dim2)
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
! generate WF on a spherical mesh
  wantp=zzero
  call timer_start(t_sic_wan_gen)
  do irloc=1,nrloc
    ir=mpi_grid_map(s_nr,dim2,loc=irloc)
    do itp=1,s_ntp
      vrc(:)=s_spx(:,itp)*s_r(ir)+wanpos(:,n)
      call s_get_wanval(n,vrc,wanval)
      wantp(itp,ir,:)=wanval(:)
    enddo
  enddo
  call mpi_grid_reduce(wantp(1,1,1),s_ntp*s_nr*nspinor,all=.true.)
! convert to spherical harmonics
  do ispn=1,nspinor
    call zgemm('T','N',lmmaxwan,s_nr,s_ntp,zone,s_ylmb,s_ntp,&
      wantp(1,1,ispn),s_ntp,zzero,s_wanlm(1,1,ispn,j),lmmaxwan)
  enddo
! compute norm
  t1=0.d0
  z1=zzero
  do ispn=1,nspinor
    do ir=1,s_nr
      do itp=1,s_ntp
        t1=t1+abs(wantp(itp,ir,ispn))**2*s_tpw(itp)*s_rw(ir)
      enddo
      z1=z1+zdotc(lmmaxwan,s_wanlm(1,ir,ispn,j),1,&
        s_wanlm(1,ir,ispn,j),1)*s_rw(ir)
    enddo
  enddo
  wanprop(wp_normtp,j)=t1
  wanprop(wp_normlm,j)=dreal(z1)
  call timer_stop(t_sic_wan_gen)
  if (irms.eq.1) then
    call timer_start(t_sic_wan_rms)
! check expansion
    call sic_wanrms(n,wantp,s_wanlm(1,1,1,j),wanprop(1,j))
    call timer_stop(t_sic_wan_rms)
  endif
! estimate the quadratic spread <r^2>-<r>^2
  x2=0.d0
  x=0.d0
  do ir=1,s_nr
    do itp=1,s_ntp
      vrc(:)=s_spx(:,itp)*s_r(ir)
      do ispn=1,nspinor
        t1=s_tpw(itp)*s_rw(ir)*abs(wantp(itp,ir,ispn))**2
        x2=x2+t1*dot_product(vrc,vrc)
        do i=1,3
          x(i)=x(i)+t1*vrc(i)
        enddo
      enddo
    enddo
  enddo
  wanprop(wp_spread,j)=x2-dot_product(x,x)
! generate potentials
  if (sic_apply(n).eq.2) then
    call timer_start(t_sic_wan_pot)
    call s_gen_pot(s_wanlm(1,1,1,j),wantp,s_wvlm(1,1,1,j),wanprop(1,j))
    call timer_stop(t_sic_wan_pot)
  endif
enddo !j
deallocate(wantp)
! check <W|\phi> matrix elements
!if (.true.) call sic_test_fvprj(fout)
! compute overlap integrals 
allocate(ovlp(sic_wantran%nwt))
ovlp=zzero
t1=0.d0
call timer_start(t_sic_wan_ovl,reset=.true.)
nwtloc=mpi_grid_map2(sic_wantran%nwt,dims=(/dim_k,dim2/))
do iloc=1,nwtloc
  i=mpi_grid_map2(sic_wantran%nwt,dims=(/dim_k,dim2/),loc=iloc)
  n=sic_wantran%iwt(1,i)
  j=sic_wantran%idxiwan(n)
  n1=sic_wantran%iwt(2,i)
  j1=sic_wantran%idxiwan(n1)
  vl(:)=sic_wantran%iwt(3:5,i)
  pos1(:)=wanpos(:,n)
  pos2(:)=wanpos(:,n1)+vl(1)*avec(:,1)+vl(2)*avec(:,2)+vl(3)*avec(:,3)
  if ((iovlp.eq.1.and.sum(abs(pos1-pos2)).lt.1d-10).or.iovlp.eq.2) then
    do ispn=1,nspinor
      ovlp(i)=ovlp(i)+s_dot_ll(pos1,pos2,s_wanlm(1,1,ispn,j),&
        s_wanlm(1,1,ispn,j1))
    enddo
  endif
  z1=ovlp(i)
  if (n.eq.n1.and.all(vl.eq.0)) then
    z1=z1-zone
  endif
  t1=max(t1,abs(z1))
enddo
call mpi_grid_reduce(ovlp(1),sic_wantran%nwt)
call mpi_grid_reduce(t1,op=op_max)
call timer_stop(t_sic_wan_ovl)
! compute energies
! note: here Hartree potential has a positive sign and XC potential 
!  has a negative sign
sic_energy_kin=0.d0
sic_epot_h=0.d0
sic_epot_xc=0.d0
do j=1,sic_wantran%nwan
  sic_energy_kin=sic_energy_kin+wanprop(wp_vsic,j)
  sic_epot_h=sic_epot_h+0.5d0*wanprop(wp_vha,j)
  sic_epot_xc=sic_epot_xc+wanprop(wp_exc,j)
enddo
sic_energy_pot=sic_epot_h+sic_epot_xc
! total energy engytot=engytot+sic_energy_tot
sic_energy_tot=sic_energy_kin-sic_energy_pot
! print some info
if (wproc) then
  write(fout,*)
  write(fout,'("   n |     normlm      normtp      spread&
    &       RMS (wan, rho, rho^{1/3})")')
  write(fout,'(80("-"))')
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    write(151,'(I4," | ",6(F10.6,2X))')n,wanprop(wp_normlm,j),&
    wanprop(wp_normtp,j),wanprop(wp_spread,j),wanprop(wp_rmswan,j),&
    wanprop(wp_rmsrho,j),wanprop(wp_rmsrho13,j)
  enddo
  write(fout,'(80("-"))')
  write(fout,'("maximum deviation from norm : ",F12.6)')t1
  write(fout,'("total quadratic spread : ",F12.6," [a.u.]^2")')&
    sum(wanprop(wp_spread,:))
  write(fout,*)
  write(fout,'("   n | ",5X,"V_n^{H}     V_n^{XC}          V_n     E_n^{XC}")')
  write(fout,'(80("-"))')
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    write(fout,'(I4," | ",4(F12.6,1X))')n,wanprop(wp_vha,j),&
      wanprop(wp_vxc,j),wanprop(wp_vsic,j),wanprop(wp_exc,j)
  enddo
  write(fout,'(80("-"))')
  write(fout,'("SIC total energy contribution      : ",G18.10,&
    &"  ! kinetic - potential ")')sic_energy_tot
  write(fout,'("SIC kinetic energy contribution    : ",G18.10,&
    &"  ! sum of V_n" )')sic_energy_kin
  write(fout,'("SIC potential energy contribution  : ",G18.10,&
    &"  ! sum of Hartree and XC terms")')sic_energy_pot
  write(fout,'("  Hartree                          : ",G18.10,&
    &"  ! sum of V_n^{H}/2")')sic_epot_h
  write(fout,'("  XC                               : ",G18.10,&
    &"  ! sum of E_n^{XC}")')sic_epot_xc 
  call flushifc(fout)
endif
deallocate(ovlp)
deallocate(wanprop)
call timer_stop(t_sic_wan)
if (wproc) then
  write(fout,*)
  write(fout,'("total time          : ",F8.3," sec.")')timer_get_value(t_sic_wan)
  write(fout,'("  generation of WFs : ",F8.3," sec.")')timer_get_value(t_sic_wan_gen)
  write(fout,'("  rms               : ",F8.3," sec.")')timer_get_value(t_sic_wan_rms)
  write(fout,'("  potential         : ",F8.3," sec.")')timer_get_value(t_sic_wan_pot)
  write(fout,'("  overlap           : ",F8.3," sec.")')timer_get_value(t_sic_wan_ovl)
  call flushifc(fout)
endif
return
end
