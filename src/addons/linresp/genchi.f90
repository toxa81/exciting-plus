subroutine genchi(iq)
use modmain
use mod_addons_q
use mod_wannier
use mod_expigqr
use mod_linresp
implicit none
! arguments
integer, intent(in) :: iq
! local variables
complex(8), allocatable :: vcwan(:,:)
! Fourier-transform of Ixc(r)=Bxc(r)/m(r)
complex(8), allocatable :: ixcft(:)
complex(8), allocatable :: fxckrnl(:,:)
character, parameter :: orb(4)=(/'s','p','d','f'/)
complex(8), allocatable :: chi0w0(:,:)

integer i,iw,j,ifxc
integer nwloc,iwloc
integer ig
character*100 qnm,qdir,fout
real(8) fxca

real(8) t1,t2,t3,t4,t5,t6,t7

call getqdir(iq,vqm(:,iq),qdir)
call getqname(vqm(:,iq),qnm)
qnm=trim(qdir)//"/"//trim(qnm)
wproc=.false.
if (mpi_grid_root((/dim_k/))) then
  wproc=.true.
  fout=trim(qnm)//"_CHI.OUT"
  open(150,file=trim(fout),form="FORMATTED",status="REPLACE")
endif

if (wproc) then
  write(150,*)
  if (lrtype.eq.0) then
    write(150,'("Calculation of charge response functions")')
  endif
  if (lrtype.eq.1) then
    write(150,'("Calculation of magnetic response functions")')  
  endif
  write(150,*)
  write(150,'("Energy mesh parameters:")')
  write(150,'("  energy interval [eV] : ", 2F9.4)')lr_w0*ha2ev,lr_w1*ha2ev
  write(150,'("  energy step     [eV] : ", F9.4)')lr_dw*ha2ev
  write(150,'("  eta             [eV] : ", F9.4)')lr_eta*ha2ev
  write(150,*)
  if (lrtype.eq.0) then
    write(150,'("Type of fxc kernel : ")')
    if (fxctype.eq.0) write(150,'("  fxc=0 (RPA)")')
    if (fxctype.eq.1) write(150,'("  fxc=-A/2 \delta_{GG''}")')
    if (fxctype.eq.2) write(150,'("  fxc=-4*Pi*A/|G+q| \delta_{GG''}")')
  endif
  write(150,*)  
  call flushifc(150)
endif

call papi_timer_start(pt_chi)
  
! for response in Wannier bais
if (wannier_chi0_chi) then
  allocate(vcwan(megqwantran%nwt,megqwantran%nwt))
! Coulomb matrix in local basis
  vcwan=zzero
  do i=1,megqwantran%nwt
    do j=1,megqwantran%nwt
      do ig=1,ngq(iq)
        vcwan(i,j)=vcwan(i,j)+dconjg(megqwan(i,ig))*megqwan(j,ig)*vhgq(ig,iq)
      enddo
    enddo
  enddo
endif !wannier_chi0_chi

allocate(fxckrnl(ngq(iq),ngq(iq)))
allocate(ixcft(ngvec))
allocate(chi0w0(ngq(iq),ngq(iq)))
if (allocated(f_response)) deallocate(f_response)
allocate(f_response(nf_response,lr_nw,nfxca))
f_response=zzero

! distribute frequency points over 1-st dimension
nwloc=mpi_grid_map(lr_nw,dim_k)

call timer_start(1,reset=.true.)
call timer_reset(2)
call timer_reset(3)
call timer_reset(4)
call timer_reset(5)
call timer_reset(6)
call timer_reset(7)
call timer_reset(8)
if (fxctype.eq.3) then
  chi0w0=zzero
  if (mpi_grid_root()) chi0w0=chi0loc(:,:,1)
  call mpi_grid_bcast(chi0w0(1,1),ngq(iq)*ngq(iq))
endif
do iwloc=1,nwloc
  iw=mpi_grid_map(lr_nw,dim_k,loc=iwloc)
! loop over fxc
  do ifxc=1,nfxca
    fxca=fxca0+(ifxc-1)*fxca1
! prepare fxc kernel
    fxckrnl=zzero
    if (lrtype.eq.0) then
      do ig=1,ngq(iq)
        if (fxctype.eq.1) then
          fxckrnl(ig,ig)=fxckrnl(ig,ig)-fxca/2.d0
        endif
        if (fxctype.eq.2) then
          fxckrnl(ig,ig)=fxckrnl(ig,ig)-fxca*vhgq(ig,iq)
        endif
      enddo
      if (fxctype.eq.3.and.ifxc.eq.2) then
        call bsfxc(iq,chi0w0,fxckrnl)
        !call bsfxc(iq,chi0loc(1,1,iwloc),fxckrnl)
      endif
    endif !lrtype.eq.0 
    call timer_start(6)
    call solve_chi(iq,lr_w(iw),chi0loc(1,1,iwloc),fxckrnl,f_response(1,iw,ifxc))
    call timer_stop(6)
    if (wannier_chi0_chi.and.ifxc.eq.1) then
      call timer_start(7)
      call solve_chi_wan(iq,lr_w(iw),vcwan,chi0wanloc(1,1,iwloc),f_response(1,iw,ifxc))
      call timer_stop(7)
    endif !wannier_chi0_chi.and.ifxc.eq.1
  enddo !ifxcloc
enddo !iwloc
t1=timer_get_value(1)
t2=timer_get_value(2)
t3=timer_get_value(3)
t4=timer_get_value(4)
t5=timer_get_value(5)
t6=timer_get_value(6)
t7=timer_get_value(7)
if (wproc) then
  write(150,*)
  write(150,'("Total time per frequency point   : ",F8.2)')t1/lr_nw
  write(150,'("  Bloch basis part (chi)         : ",F8.2)')t6/lr_nw  
  write(150,'("  Wannier basis part (chi)       : ",F8.2)')t7/lr_nw
  call flushifc(150)
endif

call mpi_grid_reduce(f_response(1,1,1),nf_response*lr_nw*nfxca,dims=(/dim_k/))
if (mpi_grid_root(dims=(/dim_k/))) then
! write response functions to .dat file
  do ifxc=1,nfxca
    call write_response_f(iq,vqm(1,iq),ifxc)
  enddo
endif

call papi_timer_stop(pt_chi)

call mpi_grid_barrier(dims=(/dim_k/))

deallocate(fxckrnl)
deallocate(ixcft)
deallocate(chi0w0)
if (wannier_chi0_chi) then
  deallocate(vcwan)
endif
deallocate(chi0loc)
if (wannier_chi0_chi) deallocate(chi0wanloc)
if (wproc) then
  write(150,*)
  write(150,'("Done.")')
  call flushifc(150)
  close(150)
endif
return
end
