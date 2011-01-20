subroutine sic_pot(fout,ene)
use modmain
use mod_sic
use mod_nrkp
use mod_addons_q
implicit none
! arguments
integer, intent(in) :: fout
complex(8), intent(out) :: ene(4,sic_wantran%nwan)
! local variables
integer ispn,i,n,j
! potential (Hartree+XC) of Wannier function charge density
real(8), allocatable :: vhxcmt(:,:,:,:,:)
real(8), allocatable :: vhxcir(:,:,:,:)
real(8) sic_epot_h,sic_epot_xc
!complex(8), allocatable :: zm1(:,:)

if (wproc) then
  write(fout,*)
  write(fout,'(80("="))')
  write(fout,'("generating potential (Hartree+XC) of Wannier functions")')
  write(fout,'(80("="))')
endif
sic_orbitals%wvmt=zzero
sic_orbitals%wvir=zzero
!-------------------!
! Hartree potential !
!-------------------!
! generate Hartree potential of Wannier functions
!   wvmt and wvir arrays are used as temporary
call sic_genvhart(sic_orbitals%wvmt,sic_orbitals%wvir)
! deallocate unnecessary arrays
deallocate(wfsvmtnrloc)
deallocate(wfsvitnrloc)
deallocate(wanncnrloc)
! restore wproc
wproc=mpi_grid_root()
if (wproc) then
  write(fout,*)
  write(fout,'("ngqmax : ",I4)')ngqmax
  write(fout,'("time for q-vectors : ",F8.3)')timer_get_value(10)
  write(fout,'("time for Hartree potential : ",F8.3)')timer_get_value(11)
  write(fout,'("maximum absolute imaginary part (mt,ir) : ",2G18.10)') &
    maxval(abs(dimag(sic_orbitals%wvmt))),maxval(abs(dimag(sic_orbitals%wvir)))
  call timestamp(fout,"done with Hartree potential")
endif
allocate(vhxcmt(lmmaxvr,nmtloc,sic_orbitals%ntr,nspinor,sic_wantran%nwan))
allocate(vhxcir(ngrloc,sic_orbitals%ntr,nspinor,sic_wantran%nwan))
do j=1,sic_wantran%nwan
  do ispn=1,nspinor
    vhxcmt(:,:,:,ispn,j)=dreal(sic_orbitals%wvmt(:,:,:,1,j))
    vhxcir(:,:,ispn,j)=dreal(sic_orbitals%wvir(:,:,1,j))
  enddo
enddo
ene=zzero
!------------------------------!
! Hartree energy <W_n|V^H|W_n> !
!------------------------------!
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  do ispn=1,nspinor
    ene(1,j)=ene(1,j)+sic_int_zdz(sic_orbitals%wanmt(1,1,1,ispn,j),&
      sic_orbitals%wanir(1,1,ispn,j),vhxcmt(1,1,1,ispn,j),vhxcir(1,1,ispn,j),&
      sic_orbitals%wanmt(1,1,1,ispn,j),sic_orbitals%wanir(1,1,ispn,j),&
      sic_orbitals%twanmtuc(1,n))
  enddo
enddo
!-------------------------!
! XC potential and energy !
!-------------------------!
call timer_start(12,reset=.true.)
call sic_genvxc(vhxcmt,vhxcir,ene)
call timer_stop(12)
if (wproc) then
  write(fout,'("time for XC potential : ",F8.3)')timer_get_value(12)
endif
! compute <W_n|V_n|W_n>
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  do ispn=1,nspinor
    ene(3,j)=ene(3,j)+sic_int_zdz(sic_orbitals%wanmt(1,1,1,ispn,j),&
      sic_orbitals%wanir(1,1,ispn,j),vhxcmt(1,1,1,ispn,j),vhxcir(1,1,ispn,j),&
      sic_orbitals%wanmt(1,1,1,ispn,j),sic_orbitals%wanir(1,1,ispn,j),&
      sic_orbitals%twanmtuc(1,n))
  enddo
enddo
! compute <W_n|V_n^{XC}|W_n>
do j=1,sic_wantran%nwan
  ene(2,j)=ene(3,j)-ene(1,j)
enddo
! note: here Hartree potential has a positive sign and XC potential 
!  has a negative sign
sic_energy_kin=0.d0
sic_epot_h=0.d0
sic_epot_xc=0.d0
do j=1,sic_wantran%nwan
  sic_energy_kin=sic_energy_kin+dreal(ene(3,j))
  sic_epot_h=sic_epot_h+0.5d0*dreal(ene(1,j))
  sic_epot_xc=sic_epot_xc+dreal(ene(4,j))
enddo
sic_energy_pot=sic_epot_h+sic_epot_xc
! total energy: engytot=engytot+sic_etot_correction
sic_energy_tot=sic_energy_kin-sic_energy_pot
if (wproc) then
  do j=1,sic_wantran%nwan
    do i=1,4
      if (abs(dimag(ene(i,j))).gt.1d-10) then
        write(fout,'("Warning : big imaginary part of energy")')
        write(fout,'(" i, j, img : ",2I4,G18.10)')i,j,dimag(ene(i,j))
      endif
    enddo
  enddo
  write(fout,*)
  write(fout,'(2X,"wann",3X,"<W_n|V_n^{H}|W_n>   <W_n|V_n^{XC}|W_n>  &
    &<W_n|V_n|W_n>     <W_n|E_n^{XC}|W_n>")')
  write(fout,'(84("-"))')
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    write(fout,'(I4,4X,4(G18.10,2X))')n,dreal(ene(1,j)),dreal(ene(2,j)),&
      dreal(ene(3,j)),dreal(ene(4,j))
  enddo
  write(fout,*)
  write(fout,'("SIC kinetic energy contribution    : ",G18.10,&
    &"  ! sum of <W_n|V_n|W_n>" )')sic_energy_kin
  write(fout,'("SIC potential energy contribution  : ",G18.10,&
    &"  ! sum of Hartree and XC terms")')sic_energy_pot
  write(fout,'("  Hartree                          : ",G18.10,&
    &"  ! sum of 1/2 <W_n|V_n^{H}|W_n>")')sic_epot_h
  write(fout,'("  XC                               : ",G18.10,&
    &"  ! sum of <W_n|E_n^{XC}|W_n>")')sic_epot_xc 
  write(fout,'("SIC total energy contribution      : ",G18.10,&
    &"  ! kinetic - potential ")')sic_energy_tot
  call flushifc(fout)
endif
sic_orbitals%wvmt=zzero
sic_orbitals%wvir=zzero
! multiply Wannier function by potential and change sign
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  do ispn=1,nspinor
    call sic_mul_zd(-zone,sic_orbitals%wanmt(1,1,1,ispn,j),&
      sic_orbitals%wanir(1,1,ispn,j),vhxcmt(1,1,1,ispn,j),vhxcir(1,1,ispn,j),&
      sic_orbitals%wvmt(1,1,1,ispn,j),sic_orbitals%wvir(1,1,ispn,j),&
      sic_orbitals%twanmtuc(1,n))    
  enddo
enddo
!!allocate(zm1(lmmaxvr,nrmtmax))
!!do n=1,nwantot
!!  call zgemm('T','N',lmmaxvr,nrmtmax,lmmaxvr,zone,dzsht,lmmaxapw,&
!!    wanmt(1,1,1,1,1,n),lmmaxvr,zzero,zm1,lmmaxvr)
!!  do lm=1,16
!!    do ir=1,nrmtmax
!!      !write(100+n,*)ir,dreal(zm1(lm,ir)),dimag(zm1(lm,ir)) !vhxcmt(lm,ir,1,1,1,n)
!!      write(200+n,*)ir,vhxcmt(lm,ir,1,1,1,n)
!!    enddo
!!    write(100+n,*)
!!    write(200+n,*)
!!  enddo
!!enddo
!!call lf_write("vxc1.dat",vhxcmt(1,1,1,1,1,1),vhxcir(1,1,1,1))
deallocate(vhxcmt,vhxcir)
return
end
