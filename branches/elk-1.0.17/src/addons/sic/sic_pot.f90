subroutine sic_pot(fout,ene)
use modmain
use mod_lf
use mod_nrkp
use mod_addons_q
implicit none
! arguments
integer, intent(in) :: fout
complex(8), intent(out) :: ene(4,nwann)
! local variables
integer nloc,ispn,i,n,lm,ir
! potential (Hartree+XC) of Wannier function charge density
real(8), allocatable :: vhxcmt(:,:,:,:,:,:)
real(8), allocatable :: vhxcir(:,:,:,:)
real(8) sic_ekin,sic_epot
complex(8), allocatable :: zm1(:,:)

if (wproc) then
  write(fout,*)
  write(fout,'("sic_pot.f90")')
  write(fout,'("generate potential (Hartree+XC) of Wannier functions")')
  write(fout,'(80("-"))')
endif
wvmt=zzero
wvir=zzero
!-------------------!
! Hartree potential !
!-------------------!
! generate Hartree potential of Wannier functions
!  use wvmt,wvir arrays as temporary
call sic_genvhart(wvmt,wvir)
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
  write(fout,'("average imaginary part (mt,ir) : ",2G18.10)') &
    sum(abs(dimag(wvmt)))/lmmaxvr/nrmtmax/natmtot/ntrloc/nwann,&
    sum(abs(dimag(wvir)))/ngrtot/ntrloc/nwann
  call timestamp(fout,"done with Hartree potential")
endif
allocate(vhxcmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor,nwannloc))
allocate(vhxcir(ngrtot,ntrloc,nspinor,nwannloc))
do ispn=1,nspinor
  do nloc=1,nwannloc
    vhxcmt(:,:,:,:,ispn,nloc)=dreal(wvmt(:,:,:,:,1,nloc))
    vhxcir(:,:,ispn,nloc)=dreal(wvir(:,:,1,nloc))
  enddo
enddo
ene=zzero
!------------------------------!
! Hartree energy <W_n|V^H|W_n> !
!------------------------------!
do nloc=1,nwannloc
  n=mpi_grid_map(nwann,dim_k,loc=nloc)
  do ispn=1,nspinor
    ene(1,n)=ene(1,n)+lf_intgr_zdz(wanmt(1,1,1,1,ispn,nloc),wanir(1,1,ispn,nloc),&
      vhxcmt(1,1,1,1,ispn,nloc),vhxcir(1,1,ispn,nloc),(/0,0,0/),wanmt(1,1,1,1,ispn,nloc),&
      wanir(1,1,ispn,nloc))
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
do nloc=1,nwannloc
  n=mpi_grid_map(nwann,dim_k,loc=nloc)
  do ispn=1,nspinor
    ene(3,n)=ene(3,n)+lf_intgr_zdz(wanmt(1,1,1,1,ispn,nloc),wanir(1,1,ispn,nloc),&
      vhxcmt(1,1,1,1,ispn,nloc),vhxcir(1,1,ispn,nloc),(/0,0,0/),wanmt(1,1,1,1,ispn,nloc),&
      wanir(1,1,ispn,nloc))
  enddo
enddo
call mpi_grid_reduce(ene(1,1),4*nwann,dims=(/dim_k/))
! compute <W_n|V_n^{XC}|W_n>
do n=1,nwann
  ene(4,n)=ene(3,n)-ene(1,n)
enddo
sic_ekin=0.d0
sic_epot=0.d0
do n=1,nwann
  sic_ekin=sic_ekin+dreal(ene(3,n))
  sic_epot=sic_epot+0.5d0*dreal(ene(1,n))+dreal(ene(2,n))
enddo
sic_etot_correction=sic_ekin-sic_epot
if (wproc) then
  do n=1,nwann
    do i=1,4
      if (abs(dimag(ene(i,n))).gt.1d-10) then
        write(fout,'("Warning : big imaginary part of energy")')
        write(fout,'(" i, n, img : ",2I4,G18.10)')i,n,dimag(ene(i,n))
      endif
    enddo
  enddo
  write(fout,*)
  write(fout,'(2X,"wann",3X,"<W_n|V_n^{H}|W_n>   <W_n|V_n^{XC}|W_n>  &
    &<W_n|V_n|W_n>     <W_n|E_n^{XC}|W_n>")')
  write(fout,'(84("-"))')
  do n=1,nwann
    write(fout,'(I4,4X,4(G18.10,2X))')n,dreal(ene(1,n)),dreal(ene(4,n)),&
      dreal(ene(3,n)),dreal(ene(2,n))
  enddo
  write(fout,*)
  write(fout,'("SIC kinetic energy contribution   : ",G18.10)')sic_ekin
  write(fout,'("SIC potential energy contribution : ",G18.10)')sic_epot
  write(fout,'("SIC total energy contribution     : ",G18.10)')sic_etot_correction
  call flushifc(fout)
endif
! multiply Wannier function by potential and change sign
do nloc=1,nwannloc
  do ispn=1,nspinor
    call lf_mult_zd(-zone,wanmt(1,1,1,1,ispn,nloc),wanir(1,1,ispn,nloc), &
      vhxcmt(1,1,1,1,ispn,nloc),vhxcir(1,1,ispn,nloc),wvmt(1,1,1,1,ispn,nloc),&
      wvir(1,1,ispn,nloc))    
  enddo
enddo
!allocate(zm1(lmmaxvr,nrmtmax))
!do n=1,nwann
!  call zgemm('T','N',lmmaxvr,nrmtmax,lmmaxvr,zone,dzsht,lmmaxapw,&
!    wanmt(1,1,1,1,1,n),lmmaxvr,zzero,zm1,lmmaxvr)
!  do lm=1,16
!    do ir=1,nrmtmax
!      !write(100+n,*)ir,dreal(zm1(lm,ir)),dimag(zm1(lm,ir)) !vhxcmt(lm,ir,1,1,1,n)
!      write(200+n,*)ir,vhxcmt(lm,ir,1,1,1,n)
!    enddo
!    write(100+n,*)
!    write(200+n,*)
!  enddo
!enddo
deallocate(vhxcmt,vhxcir)
return
end