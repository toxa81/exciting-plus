#ifdef _HDF5_
subroutine genchi(ivq0m)
use modmain
implicit none
integer, intent(in) :: ivq0m(3)

complex(8), allocatable :: chi0(:,:)
complex(8), allocatable :: chi0wan(:,:)
complex(8), allocatable :: vcwan(:,:)
! Fourier-transform of Ixc(r)=Bxc(r)/m(r)
complex(8), allocatable :: ixcft(:)
complex(8), allocatable :: krnl(:,:)
complex(8), allocatable :: krnl_scr(:,:)
real(8), allocatable :: vcgq(:)

! G+q vector in Cartesian coordinates
real(8) vgq0c(3)
! length of G+q vector
real(8) gq0

real(8) fxca
complex(8) zt1
real(8) dvec(3),pos1(3),pos2(3),vtrc(3)
integer ias1,ias2
real(8) d

logical exist

integer ig,i,j,ig1,ig2,idx0,bs,n1,n2,it1,it2,n3,n4,n
integer i1,i2,ifxc,ifxc1,ifxc2
integer iv(3)
character*100 fchi,fchi0,fme,qnm,path,fw
character*3 c3
integer nwstep,iwstep,iw,nwloc,iwloc

! we need Bcx and magnetization from STATE.OUT
if (lrtype.eq.1) call readstate

call qname(ivq0m,qnm)
qnm="./"//trim(qnm)//"/"//trim(qnm)
wproc=.false.

fchi0=trim(qnm)//"_chi0.hdf5"

allocate(chi0(ngvecme,ngvecme))
allocate(krnl(ngvecme,ngvecme))
!if (screened_w) allocate(krnl_scr(ngvecchi,ngvecchi))
allocate(ixcft(ngvec))
if (allocated(f_response)) deallocate(f_response)
allocate(f_response(nf_response,nepts,nfxca))
f_response=zzero

! setup sqrt(4Pi)/|G+q| array
allocate(vcgq(ngvecme))
do ig=1,ngvecme
! generate G+q vectors  
  vgq0c(:)=vgc(:,ig+gvecme1-1)+vq0rc(:)
  gq0=sqrt(vgq0c(1)**2+vgq0c(2)**2+vgq0c(3)**2)
  vcgq(ig)=sqrt(fourpi)/gq0
enddo !ig


! for response in Wannier basis
if (wannier_chi0_chi) then
  allocate(chi0wan(nmegqwan,nmegqwan))
  allocate(vcwan(nmegqwan,nmegqwan))
! Coulomb matrix in local basis
  vcwan=zzero
  do i=1,nmegqwan
    do j=1,nmegqwan
      do ig=1,ngvecme
        vcwan(i,j)=vcwan(i,j)+dconjg(megqwan(i,ig))*megqwan(j,ig)*vcgq(ig)**2
      enddo
    enddo
  enddo
endif !wannier_chi0_chi

! distribute energy points between 1-st dimension
!i=0
!nwstep=mpi_grid_map(nepts,dim_k,x=i)
nwloc=mpi_grid_map(nepts,dim_k)
! distribute nfxca between 2-nd dimension 
bs=mpi_grid_map(nfxca,dim_b,offs=idx0)
ifxc1=idx0+1
ifxc2=idx0+bs
! main loop over energy points 
do iwloc=1,nwloc
  iw=mpi_grid_map(nepts,dim_k,loc=iwloc)
  if (mpi_grid_root(dims=(/dim_b/))) then
    write(path,'("/iw/",I8.8)')iw
    call read_real8_array(chi0,3,(/2,ngvecme,ngvecme/),trim(fchi0),&
      trim(path),'chi0')
  endif
  call mpi_grid_bcast(chi0(1,1),ngvecme*ngvecme,dims=(/dim_b/))
  if (crpa.and.mpi_grid_root(dims=(/dim_b/))) then
    call genwu(ngvecme,iw,chi0,vcgq,qnm)
  endif
  if (.not.crpa) then
    do ifxc=ifxc1,ifxc2
      fxca=fxca0+(ifxc-1)*fxca1
! prepare fxc kernel
      krnl=zzero
      if (lrtype.eq.0) then
        do ig=1,ngvecme
          if (fxctype.eq.1) then
            krnl(ig,ig)=krnl(ig,ig)-fxca/2.d0
          endif
          if (fxctype.eq.2) then
            krnl(ig,ig)=krnl(ig,ig)-fxca*vcgq(ig)**2
          endif
        enddo
      endif !lrtype.eq.0 
      call timer_start(6)
      call solve_chi(vcgq,lr_w(iw),chi0,krnl,krnl_scr,f_response(1,iw,ifxc))
      call timer_stop(6)
      if (wannier_chi0_chi.and.ifxc.eq.1) then
        call timer_start(7)
        call solve_chi_wan(vcgq,lr_w(iw),vcwan,chi0wan,f_response(1,iw,ifxc))
        call timer_stop(7)
      endif
    enddo !ifxc
  endif ! .not.crpa
enddo
call mpi_grid_reduce(f_response(1,1,1),nf_response*nepts*nfxca,dims=(/dim_k,dim_b/))
! write response functions to .dat file
if (mpi_grid_root(dims=(/dim_k,dim_b/))) then
  do ifxc=1,nfxca
    call write_chi(ivq0m,ifxc)
  enddo
endif
deallocate(chi0)
deallocate(krnl)
deallocate(ixcft)
deallocate(vcgq)
if (wannier_chi0_chi) then
  deallocate(chi0wan)
  deallocate(vcwan)
endif
return
end  
#endif
