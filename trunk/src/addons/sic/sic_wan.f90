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
integer n,ispn,vl(3),n1,i,j,j1,itp,ir,nwtloc,iloc,ntp,nr
real(8) t1,t2,vrc(3),x(3)
real(8) sic_epot_h,sic_epot_xc
complex(8) z1,zt1,wanval(nspinor)
real(8), allocatable :: ene(:,:)
real(8), allocatable :: xmt(:,:,:,:)
real(8), allocatable :: xir(:,:,:)
real(8), allocatable :: spread(:,:)
real(8), allocatable :: tp(:,:)
complex(8), allocatable :: ovlp(:)
complex(8), allocatable :: wantp(:,:,:)
real(8), allocatable :: wanrms(:)

allocate(wantp(s_ntp,s_nr,nspinor))
allocate(ene(4,sic_wantran%nwan))
allocate(wanrms(sic_wantran%nwan))

if (wproc) then
  write(fout,*)
  write(fout,'(80("="))')
  write(fout,'("generating Wannier functions on a spherical grid")')
  write(fout,'(80("="))')
endif
call timer_start(3,reset=.true.)
ntp=250
nr=300
allocate(tp(2,ntp))
call sphcover(ntp,tp)
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  do ir=1,s_nr
    do itp=1,s_ntp
      vrc(:)=s_spx(:,itp)*s_r(ir)
      call s_get_wanval(n,vrc,wanval)
      wantp(itp,ir,:)=wanval(:)
    enddo
  enddo
  call mpi_grid_reduce(wantp(1,1,1),s_ntp*s_nr*nspinor,dims=(/dim_k/),all=.true.)
! convert to spherical harmonics
  do ispn=1,nspinor
    call zgemm('T','N',lmmaxwan,s_nr,s_ntp,zone,s_ylmb,s_ntp,wantp(1,1,ispn),&
      s_ntp,zzero,s_wanlm(1,1,ispn,j),lmmaxwan)
  enddo
! check expansion
  t1=0.d0
  do ir=1,nr
    do itp=1,ntp
      vrc(:)=(/sin(tp(1,itp))*cos(tp(2,itp)),&
               sin(tp(1,itp))*sin(tp(2,itp)),&
               cos(tp(1,itp))/)*ir*sic_wan_cutoff/nr/2.d0
      call s_get_wanval(n,vrc,wanval)
      do ispn=1,nspinor
        z1=wanval(ispn)-s_func_val(vrc,s_wanlm(1,1,ispn,j))
        t1=t1+abs(z1)**2
      enddo
    enddo
  enddo
  wanrms(j)=sqrt(t1/nr/ntp)
! generate potentials
  if (sic_apply(n).eq.2) then
    call s_gen_pot(s_wanlm(1,1,1,j),wantp,s_wvlm(1,1,1,j),ene(1,j),ene(2,j),&
      ene(3,j),ene(4,j))
  else
    s_wvlm(:,:,:,j)=zzero
    ene(:,j)=0.d0
  endif
enddo !j
deallocate(wantp)
deallocate(tp)

allocate(ovlp(sic_wantran%nwt))
ovlp=zzero
t1=0.d0
! compute overlap integrals 
nwtloc=mpi_grid_map(sic_wantran%nwt,dim_k)
do iloc=1,nwtloc
  i=mpi_grid_map(sic_wantran%nwt,dim_k,loc=iloc)
  n=sic_wantran%iwt(1,i)
  j=sic_wantran%idxiwan(n)
  n1=sic_wantran%iwt(2,i)
  j1=sic_wantran%idxiwan(n1)
  vl(:)=sic_wantran%iwt(3:5,i)
  do ispn=1,nspinor
    ovlp(i)=ovlp(i)+s_dot_ll(n,n1,vl,s_wanlm(1,1,ispn,j),s_wanlm(1,1,ispn,j1))
  enddo
  z1=ovlp(i)
  if (n.eq.n1.and.all(vl.eq.0)) then
    z1=z1-zone
  endif
  t1=max(t1,abs(z1))
enddo
call mpi_grid_reduce(ovlp(1),sic_wantran%nwt,dims=(/dim_k/))
call mpi_grid_reduce(t1,dims=(/dim_k/),op=op_max)
! compute energies
! note: here Hartree potential has a positive sign and XC potential 
!  has a negative sign
sic_energy_kin=0.d0
sic_epot_h=0.d0
sic_epot_xc=0.d0
do j=1,sic_wantran%nwan
  sic_energy_kin=sic_energy_kin+ene(3,j)
  sic_epot_h=sic_epot_h+0.5d0*ene(1,j)
  sic_epot_xc=sic_epot_xc+ene(4,j)
enddo
sic_energy_pot=sic_epot_h+sic_epot_xc
! total energy: engytot=engytot+sic_etot_correction
sic_energy_tot=sic_energy_kin-sic_energy_pot
! print some info
if (wproc) then
  write(fout,*)
  write(fout,'(" n   |   norm     RMS     spread")')
  write(fout,'(60("-"))')
  do i=1,sic_wantran%nwan
    n=sic_wantran%iwan(i)
    j=sic_wantran%iwtidx(n,n,0,0,0)
    write(151,'(I4," | ",2X,G18.10,2X,G18.10)')n,dreal(ovlp(j)),wanrms(i)
  enddo
  write(fout,'(60("-"))')
  write(fout,'("maximum deviation from norm : ",F12.6)')t1
  write(fout,*)
  write(fout,'(2X,"n  | ",3X,"V_n^{H}  V_n^{XC} V_n  E_n^{XC}")')
  write(fout,'(84("-"))')
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    write(fout,'(I4," | ",4X,4(G18.10,2X))')n,ene(:,j)
  enddo
  write(fout,'(84("-"))')
  write(fout,*)
  write(fout,'("SIC kinetic energy contribution    : ",G18.10,&
    &"  ! sum of V_n" )')sic_energy_kin
  write(fout,'("SIC potential energy contribution  : ",G18.10,&
    &"  ! sum of Hartree and XC terms")')sic_energy_pot
  write(fout,'("  Hartree                          : ",G18.10,&
    &"  ! sum of V_n^{H}/2")')sic_epot_h
  write(fout,'("  XC                               : ",G18.10,&
    &"  ! sum of E_n^{XC}")')sic_epot_xc 
  write(fout,'("SIC total energy contribution      : ",G18.10,&
    &"  ! kinetic - potential ")')sic_energy_tot
  call flushifc(fout)
endif
deallocate(ovlp)
deallocate(wanrms)

!allocate(xmt(lmmaxvr,nmtloc,sic_orbitals%ntr,4))
!allocate(xir(ngrloc,sic_orbitals%ntr,4))
!call sic_gen_r(xmt,xir)
!allocate(spread(sic_wantran%nwan,6))
!spread=0.d0
!do j=1,sic_wantran%nwan
!  n=sic_wantran%iwan(j)
!  z2=zzero
!  do i=1,4
!    do ispn=1,nspinor
!      z2(i)=z2(i)+sic_int_zdz(sic_orbitals%wanmt(1,1,1,ispn,j),&
!        sic_orbitals%wanir(1,1,ispn,j),xmt(1,1,1,i),xir(1,1,i),&
!        sic_orbitals%wanmt(1,1,1,ispn,j),sic_orbitals%wanir(1,1,ispn,j),&
!        sic_orbitals%twanmtuc(1,n))
!    enddo
!  enddo
!  do i=1,3
!    vrc(i)=dreal(z2(i))
!  enddo
!  spread(j,1)=dreal(z2(4))-dot_product(vrc(:),vrc(:))
!  spread(j,2)=dreal(z2(4))
!  spread(j,3)=dot_product(vrc(:),vrc(:))
!  spread(j,4:6)=vrc(:)
!enddo
!if (wproc) then
!  write(fout,*)
!  write(fout,'("quadratic spreads",3X,"<r^2> - <r>^2  [a.u.]^2&
!  &         (<r^2>  <r>^2  <x> <y> <z>) ")')
!  write(fout,'(90("-"))')
!  do j=1,sic_wantran%nwan
!    n=sic_wantran%iwan(j)
!    write(151,'("  n : ",I4,8X,G18.10,3X,"(",5F10.5,")")')n,spread(j,:)
!  enddo
!  write(fout,'(90("-"))')
!  write(fout,'("total spread : ",F12.6)')sum(spread)
!  call timestamp(fout,"done with Wannier functions")
!  call flushifc(151)  
!endif
!deallocate(xmt,xir,spread)
call timer_stop(3)
if (wproc) then
  write(fout,*)
  write(fout,'("done in ",F8.3," sec.")')timer_get_value(3)
  call flushifc(fout)
endif
deallocate(ene)
return
end
