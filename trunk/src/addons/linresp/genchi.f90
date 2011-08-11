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
complex(8), allocatable :: krnl(:,:)
character, parameter :: orb(4)=(/'s','p','d','f'/)

integer i,iw,j,ifxc
integer nfxcloc,ifxcloc,nwloc,iwloc
integer ig
character*100 qnm,qdir,fout
real(8) fxca

real(8) t1,t2,t3,t4,t5,t6,t7

call getqdir(iq,vqm(:,iq),qdir)
call getqname(vqm(:,iq),qnm)
qnm=trim(qdir)//"/"//trim(qnm)
wproc=.false.
if (mpi_grid_root((/dim_k,dim_b/))) then
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
  write(150,'("  energy interval [eV] : ", 2F9.4)')lr_w0,lr_w1
  write(150,'("  energy step     [eV] : ", F9.4)')lr_dw
  write(150,'("  eta             [eV] : ", F9.4)')lr_eta
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
      do ig=1,ngvecme
        vcwan(i,j)=vcwan(i,j)+dconjg(megqwan(i,ig))*megqwan(j,ig)*vhgq(ig,iq)
      enddo
    enddo
  enddo
endif !wannier_chi0_chi

allocate(krnl(ngvecme,ngvecme))
allocate(ixcft(ngvec))
if (allocated(f_response)) deallocate(f_response)
allocate(f_response(nf_response,lr_nw,nfxca))
f_response=zzero

! distribute nfxca between 2-nd dimension 
nfxcloc=mpi_grid_map(nfxca,dim_b)
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
do iwloc=1,nwloc
  iw=mpi_grid_map(lr_nw,dim_k,loc=iwloc)
! broadcast chi0
  call mpi_grid_bcast(chi0loc(1,1,iwloc),ngvecme*ngvecme,dims=(/dim_b/))
! loop over fxc
  do ifxcloc=1,nfxcloc
    ifxc=mpi_grid_map(nfxca,dim_b,loc=ifxcloc)
    fxca=fxca0+(ifxc-1)*fxca1
! prepare fxc kernel
    krnl=zzero
    if (lrtype.eq.0) then
      do ig=1,ngvecme
        if (fxctype.eq.1) then
          krnl(ig,ig)=krnl(ig,ig)-fxca/2.d0
        endif
        if (fxctype.eq.2) then
          krnl(ig,ig)=krnl(ig,ig)-fxca*vhgq(ig,iq)
        endif
      enddo
      if (fxctype.eq.3) then
        call bsfxc(iq,chi0loc(1,1,1),krnl)
      endif
    endif !lrtype.eq.0 
    call timer_start(6)
    call solve_chi(iq,lr_w(iw),chi0loc(1,1,iwloc),krnl,f_response(1,iw,ifxc))
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

call mpi_grid_reduce(f_response(1,1,1),nf_response*lr_nw*nfxca,dims=(/dim_k,dim_b/))
if (mpi_grid_root(dims=(/dim_k,dim_b/))) then
! write response functions to .dat file
  do ifxc=1,nfxca
    call write_chi(iq,vqm(1,iq),ifxc)
  enddo
endif

call papi_timer_stop(pt_chi)

call mpi_grid_barrier(dims=(/dim_k,dim_b/))

deallocate(krnl)
deallocate(ixcft)
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
