#ifdef _HDF5_
subroutine genchi0(ivq0m)
use modmain
use hdf5
implicit none
! arguments
integer, intent(in) :: ivq0m(3)
! local variables
complex(8), allocatable :: chi0(:,:)
complex(8), allocatable :: chi0w(:,:,:)
complex(8), allocatable :: chi0wan(:,:)
complex(8), allocatable :: chi0wan_k(:,:,:)
complex(8), allocatable :: vcwan(:,:)
! Fourier-transform of Ixc(r)=Bxc(r)/m(r)
complex(8), allocatable :: ixcft(:)
complex(8), allocatable :: krnl(:,:)
complex(8), allocatable :: krnl_scr(:,:)
complex(8), allocatable :: mexp(:,:,:)
integer, external :: hash


integer i,iw,i1,i2,ikloc,n,j,ifxc
integer ist1,ist2,nfxcloc,ifxcloc,nwloc,jwloc,iwloc,tag
integer it1(3),it2(3),it(3)
integer ig
character*100 qnm,fout,fchi0,fstat,path
character*8 c8
integer ie1,n1,n2,ik
real(8) fxca

real(8), allocatable :: vcgq(:)
real(8) vgq0c(3)
real(8) gq0

real(8) vtrc(3)
real(8) t1,t2,t3,t4,t5,t6,t7,t8

call qname(ivq0m,qnm)
qnm="./"//trim(qnm)//"/"//trim(qnm)
wproc=.false.
if (mpi_grid_root((/dim_k,dim_b/))) then
  wproc=.true.
  fout=trim(qnm)//"_LR.OUT"
  open(150,file=trim(fout),form='formatted',status='replace')
  fstat=trim(qnm)//"_chi0_stat.txt"
endif

! we need Bcx and magnetization from STATE.OUT
if (lrtype.eq.1) call readstate

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
  write(150,'("  maximum energy [eV] : ", F9.4)')maxomega
  write(150,'("  energy step    [eV] : ", F9.4)')domega
  write(150,'("  eta            [eV] : ", F9.4)')lr_eta
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
  
! setup sqrt(4Pi)/|G+q| array
allocate(vcgq(ngvecme))
do ig=1,ngvecme
! generate G+q vectors  
  vgq0c(:)=vgc(:,ig+gvecme1-1)+vq0rc(:)
  gq0=sqrt(vgq0c(1)**2+vgq0c(2)**2+vgq0c(3)**2)
  vcgq(ig)=sqrt(fourpi)/gq0
enddo !ig

! read matrix elements
if (write_megq_file) then
  call timer_start(1,reset=.true.)
  if (wproc) then
    write(150,*)
    write(150,'("Reading matrix elements")')
    call flushifc(150)
  endif
  call readmegqblh(qnm)
  if (wannier_chi0_chi) call readmegqwan(qnm)
  call timer_stop(1)
  if (wproc) then
    write(150,'("Done in ",F8.2," seconds")')timer_get_value(1)
    write(150,*)
    write(150,'("matrix elements were calculated for: ")')
    write(150,'("  G-shells  : ",I4," to ", I4)')gshme1,gshme2
    write(150,'("  G-vectors : ",I4," to ", I4)')gvecme1,gvecme2
    call flushifc(150)
  endif
endif
allocate(megqblh2(nmegqblhlocmax,ngvecme))

! for response in Wannier bais
if (wannier_chi0_chi) then
  if (wproc) then
    write(150,*)
    write(150,'("Wannier AFM : ",L1)')megqwan_afm
  endif
  allocate(chi0wan(nmegqwan,nmegqwan))
  allocate(chi0wan_k(nmegqwan,nmegqwan,nkptnrloc))
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
  allocate(mexp(nmegqwan,nmegqwan,nkptnrloc))
  do i1=1,nmegqwan
    do i2=1,nmegqwan
      it1(:)=imegqwan(3:5,i1)
      it2(:)=imegqwan(3:5,i2)
      it(:)=it1(:)-it2(:)
      vtrc(:)=avec(:,1)*it(1)+avec(:,2)*it(2)+avec(:,3)*it(3)
      do ikloc=1,nkptnrloc
        ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
! phase e^{i(k+q)T}
        mexp(i1,i2,ikloc)=exp(dcmplx(0.d0,dot_product(vkcnr(:,ik)+vq0rc(:),vtrc(:))))   
      enddo
    enddo
  enddo
! arrangement for zgemm  
  allocate(wann_cc(nmegqblhwanmax,nmegqwan,nkptnrloc))
  allocate(wann_cc2(nmegqblhwanmax,nmegqwan))
! arrangement for zgerc  
!  allocate(wann_cc(nmegqwan,nmegqblhwanmax,nkptnrloc))
  wann_cc=zzero
  do ikloc=1,nkptnrloc
    do i1=1,nmegqblhwan(ikloc)
      i=imegqblhwan(i1,ikloc)
      ist1=bmegqblh(1,i,ikloc)
      ist2=bmegqblh(2,i,ikloc)
      do n=1,nmegqwan
        n1=imegqwan(1,n)
        n2=imegqwan(2,n)
        ! for zgemm
        wann_cc(i1,n,ikloc)=wann_c(n1,ist1,ikloc)*dconjg(wann_c(n2,ist2,ikloc+nkptnrloc))
        ! for zgerc
        !wann_cc(n,i1,ikloc)=wann_c(n1,ist1,ikloc)*dconjg(wann_c(n2,ist2,ikloc+nkptnrloc))
      enddo
    enddo !i1
  enddo !ikloc
endif !wannier_chi0_chi

ie1=1

allocate(chi0(ngvecme,ngvecme))
allocate(krnl(ngvecme,ngvecme))
if (screened_w.or.crpa) allocate(krnl_scr(ngvecme,ngvecme))
allocate(ixcft(ngvec))
if (allocated(f_response)) deallocate(f_response)
allocate(f_response(nf_response,nepts,nfxca))
f_response=zzero

fchi0=trim(qnm)//"_chi0.hdf5"
if (mpi_grid_root((/dim_k,dim_b/)).and.write_chi0_file) then
  call hdf5_create_file(trim(fchi0))
  call hdf5_create_group(trim(fchi0),'/','iw')
  do i=1,nepts
    write(c8,'(I8.8)')i
    call hdf5_create_group(trim(fchi0),'/iw',c8)
  enddo
endif
call mpi_grid_barrier(dims=(/dim_k,dim_b/))

! distribute nfxca between 2-nd dimension 
nfxcloc=mpi_grid_map(nfxca,dim_b)
nwloc=mpi_grid_map(nepts,dim_k)

if (crpa) then
  allocate(chi0w(ngvecme,ngvecme,nwloc))
endif

if (wproc) then
  write(150,*)
  write(150,'("fisrt energy point : ",I4)')ie1
  call flushifc(150)
endif
call timer_start(1,reset=.true.)
call timer_reset(2)
call timer_reset(3)
call timer_reset(4)
call timer_reset(5)
call timer_reset(6)
call timer_reset(7)
call timer_reset(8)
! loop over energy points
do iw=ie1,nepts
  chi0=zzero
  if (wannier_chi0_chi) chi0wan_k=zzero
! sum over fraction of k-points
  call timer_start(2)
  do ikloc=1,nkptnrloc
    if (nmegqblhloc(1,ikloc).gt.0) then
! for each k-point : sum over interband transitions
      call sumchi0(ikloc,lr_w(iw),chi0)
    endif
  enddo
! sum over k-points and band transitions
  call mpi_grid_reduce(chi0(1,1),ngvecme*ngvecme,dims=(/dim_k,dim_b/))
  chi0=chi0/nkptnr/omega
  call timer_stop(2)
! for response in Wannier basis
  if (wannier_chi0_chi) then
    call timer_start(3)
    do ikloc=1,nkptnrloc
      if (nmegqblhwan(ikloc).gt.0) then
        call sumchi0wan_k(ikloc,lr_w(iw),chi0wan_k(1,1,ikloc))
      endif
    enddo !ikloc
    !call mpi_grid_reduce(chi0wan_k(1,1,1),nmegqwan*nmegqwan*nkptnrloc,&
    !  dims=(/dim_b/))
    call timer_stop(3)
    call timer_start(4)
! compute ch0 matrix in Wannier basis
    call genchi0wan(mexp,chi0wan_k,chi0wan)
    if (megqwan_afm) chi0wan(:,:)=chi0wan(:,:)*2.d0    
    call timer_stop(4)
  endif
  
  if (write_chi0_file) then
    call timer_start(8)
    if (mpi_grid_root((/dim_k,dim_b/))) then
      write(path,'("/iw/",I8.8)')iw
      call hdf5_write(trim(fchi0),trim(path),'w',lr_w(iw))
      call hdf5_write(trim(fchi0),trim(path),'chi0',chi0(1,1),&
        dims=(/ngvecme,ngvecme/))
      if (wannier_chi0_chi) then
        call hdf5_write(trim(fchi0),trim(path),'chi0wan',chi0wan(1,1), &
          dims=(/nmegqwan,nmegqwan/))
      endif
    endif
    call timer_stop(8)
  else
    if (crpa.and.mpi_grid_root(dims=(/dim_b/))) then
      jwloc=mpi_grid_map(nepts,dim_k,glob=iw,x=j)
      if (mpi_grid_x(dim_k).eq.0) then
        if (j.eq.0) then
          chi0w(:,:,jwloc)=chi0(:,:)  
        else
          tag=iw
          call mpi_grid_send(chi0(1,1),ngvecme*ngvecme,(/dim_k/),(/j/),tag)
        endif
      endif
      if (mpi_grid_x(dim_k).eq.j.and.j.ne.0) then
        tag=iw
        call mpi_grid_recieve(chi0w(1,1,jwloc),ngvecme*ngvecme,(/dim_k/),&
          (/0/),tag)
      endif
    endif
! compute response functions
    if (mpi_grid_root(dims=(/dim_k/)).and..not.crpa) then
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
      enddo
    endif !(mpi_grid_root(dims=(/dim_k/)).and..not.crpa)
  endif !write_chi0_file
  if (wproc) then
    open(160,file=trim(fstat),status='replace',form='formatted')
    write(160,'(I8)')iw
    close(160)
  endif
enddo !iw
call timer_stop(1)
do iwloc=1,nwloc
  iw=mpi_grid_map(nepts,dim_k,loc=iwloc)
  call timer_start(5)
  call genwu(iw,chi0w(1,1,iwloc),vcgq,qnm,krnl_scr)
  call timer_stop(5)
enddo
t1=timer_get_value(1)
t2=timer_get_value(2)
t3=timer_get_value(3)
t4=timer_get_value(4)
t5=timer_get_value(5)
t6=timer_get_value(6)
t7=timer_get_value(7)
t8=timer_get_value(8)
if (wproc) then
  write(150,*)
  write(150,'("Total time per frequency point   : ",F8.2)')t1/nepts
  write(150,'("  Bloch basis part (chi0)        : ",F8.2)')t2/nepts
  write(150,'("  Bloch basis part (chi)         : ",F8.2)')t6/nepts  
  write(150,'("  Wannier basis part (chi0wan_k) : ",F8.2)')t3/nepts
  write(150,'("  Wannier basis part (chi0wan)   : ",F8.2)')t4/nepts 
  write(150,'("  Wannier basis part (crpa)      : ",F8.2)')t5/nwloc   
  write(150,'("  Wannier basis part (chi)       : ",F8.2)')t7/nepts
  write(150,'("  Write chi0                     : ",F8.2)')t8/nepts  
  call flushifc(150)
endif

if (.not.write_chi0_file) then
  if (mpi_grid_root(dims=(/dim_k/))) then
    call mpi_grid_reduce(f_response(1,1,1),nf_response*nepts*nfxca,dims=(/dim_b/))
! write response functions to .dat file
    if (mpi_grid_root(dims=(/dim_b/))) then
      do ifxc=1,nfxca
        call write_chi(ivq0m,ifxc)
      enddo
    endif
  endif
endif

call mpi_grid_barrier(dims=(/dim_k,dim_b/))

deallocate(chi0)
deallocate(krnl)
deallocate(ixcft)
if (wannier_chi0_chi) then
  deallocate(chi0wan)
  deallocate(chi0wan_k)
  deallocate(vcwan)
  deallocate(mexp)
  deallocate(wann_cc)
  deallocate(wann_cc2)
endif
deallocate(megqblh2)
deallocate(vcgq)
if (crpa) then
  deallocate(chi0w)
  deallocate(krnl_scr)
endif

if (write_chi0_file) then
  call genchi(ivq0m)
endif

30 continue
if (wproc) then
  write(150,*)
  write(150,'("Done.")')
  call flushifc(150)
endif


return
end

#endif
