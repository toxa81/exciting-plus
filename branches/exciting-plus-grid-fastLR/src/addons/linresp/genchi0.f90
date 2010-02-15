#ifdef _HDF5_
subroutine genchi0(ivq0m)
use modmain
use hdf5
implicit none
! arguments
integer, intent(in) :: ivq0m(3)
! local variables
complex(8), allocatable :: chi0(:,:)
complex(8), allocatable :: chi0wan(:,:,:)
complex(8), allocatable :: chi0wan_k(:,:,:)
complex(8), allocatable :: vcwan(:,:)
complex(8) :: chi0_GqGq_wan_full
! Fourier-transform of Ixc(r)=Bxc(r)/m(r)
complex(8), allocatable :: ixcft(:)
complex(8), allocatable :: krnl(:,:)
complex(8), allocatable :: krnl_scr(:,:)
integer, external :: hash


integer i,ik,ie,i1,i2,i3,ikloc,it1,n,j,bs,ifxc1,ifxc2,ifxc,idx0
integer ig,igq0
character*100 path,qnm,fout,fchi0,fu
logical exist
integer ie1,n1,n2,jk
integer iv(3)
real(8) fxca

real(8), allocatable :: vcgq(:)
real(8) vgq0c(3)
real(8) gq0

real(8) dvec(3),pos1(3),pos2(3),vtrc(3)
integer ias1,ias2
real(8) d


logical lafm

! non-zero matrix elements in Wannier basis
integer nmegqwan2
integer, allocatable :: imegqwan2(:,:)
complex(8), allocatable :: megqwan2(:,:,:)

call qname(ivq0m,qnm)
qnm="./"//trim(qnm)//"/"//trim(qnm)
wproc=.false.
if (mpi_grid_root((/dim_k,dim_b/))) then
  wproc=.true.
  fout=trim(qnm)//"_LR.OUT"
  open(150,file=trim(fout),form='formatted',status='replace')
endif

if (crpa) then
  if (mpi_grid_root((/dim_k,dim2/))) then
    fu=trim(qnm)//"_U"
    inquire(file=trim(fu),exist=exist)
  endif
  call mpi_grid_bcast(exist,dims=(/dim_k,dim2/))
  if (exist) goto 30
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

!read matrix elements
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
  lafm=.false.
!  if ((all(iwann(3,:).eq.1).or.all(iwann(3,:).eq.2)).and.spinpol) then
!    lafm=.true.
!  endif
  if (wproc) then
    write(150,*)
    write(150,'("AFM case : ",L1)')lafm
  endif
  allocate(chi0wan(nmegqwan,nmegqwan,ntrchi0wan))
  allocate(chi0wan_k(nmegqwan,nmegqwan,nkptnrloc))
  allocate(megqwan2(nmegqwan,ntrmegqwan,ngvecme))
  megqwan2=zzero
  !if (wproc) then
  !  write(150,'("minimal matrix element : ",F12.6)')megqwan_cutoff
  !endif
  do it1=1,ntrmegqwan
    vtrc(:)=avec(:,1)*itrmegqwan(1,it1)+&
            avec(:,2)*itrmegqwan(2,it1)+&
            avec(:,3)*itrmegqwan(3,it1)
    do n=1,nmegqwan
      n1=bmegqwan(1,n)
      n2=bmegqwan(2,n)
      ias1=iwann(1,n1)
      ias2=iwann(1,n2)
      pos1(:)=atposc(:,ias2ia(ias1),ias2is(ias1))
      pos2(:)=atposc(:,ias2ia(ias2),ias2is(ias2))
      dvec(:)=pos2(:)+vtrc(:)-pos1(:)
      d=sqrt(dvec(1)**2+dvec(2)**2+dvec(3)**2)
      if (d.lt.megqwan_maxdist) megqwan2(n,it1,:)=megqwan(n,it1,:)
      !do ig=1,ngvecme
      !  if (abs(megqwan(n,it1,ig)).lt.megqwan_cutoff) megqwan(n,it1,ig)=zzero
      !enddo
    enddo
  enddo
!  inquire(file='mewf.in',exist=exist)
!  if (exist) then
!    open(70,file='mewf.in',form='formatted',status='old')
!    read(70,*)itrans_m
!    if (itrans_m.eq.1) then
!      read(70,*)d1
!      if (wproc) then
!        write(150,'("minimal matrix element : ",F12.6)')d1
!      endif
!      do it1=1,ntrmegqwan
!        do n=1,nmegqwan
!          do ig=1,ngvecme
!            if (abs(megqwan(n,it1,ig)).lt.d1) megqwan(n,it1,ig)=zzero
!          enddo
!        enddo
!      enddo
!    endif
!    if (itrans_m.eq.2) then
!      read(70,*)ntrans
!      allocate(itrans(5,ntrans))
!      do i=1,ntrans
!        read(70,*)itrans(1:5,i)
!      enddo
!      allocate(megqwan_t(nmegqwan,ntrmegqwan,ngvecme))
!      do it1=1,ntrmegqwan
!        do n=1,nmegqwan
!          n1=bmegqwan(1,n)
!          n2=bmegqwan(2,n)
!          do i=1,ntrans
!            if (n1.eq.itrans(1,i).and.&
!                n2.eq.itrans(2,i).and.&
!                itrmegqwan(1,it1).eq.itrans(3,i).and.&
!                itrmegqwan(2,it1).eq.itrans(4,i).and.&
!                itrmegqwan(3,it1).eq.itrans(5,i)) then
!              megqwan_t(n,it1,:)=megqwan(n,it1,:)
!            endif
!            if (n2.eq.itrans(1,i).and.&
!                n1.eq.itrans(2,i).and.&
!                itrmegqwan(1,it1).eq.-itrans(3,i).and.&
!                itrmegqwan(2,it1).eq.-itrans(4,i).and.&
!                itrmegqwan(3,it1).eq.-itrans(5,i)) then
!              megqwan_t(n,it1,:)=megqwan(n,it1,:)
!            endif
!          enddo
!        enddo
!      enddo
!      megqwan=megqwan_t
!      deallocate(megqwan_t)
!      deallocate(itrans)
!    endif
!  endif
  if (wproc) then
    write(150,*)
    write(150,'("Matrix elements in WF basis : ")')
    write(150,*)
    do it1=1,ntrmegqwan
      vtrc(:)=avec(:,1)*itrmegqwan(1,it1)+&
              avec(:,2)*itrmegqwan(2,it1)+&
              avec(:,3)*itrmegqwan(3,it1)    
      if (sum(abs(megqwan2(:,it1,:))).gt.1d-16) then
        write(150,'("translation : ",3I4)')itrmegqwan(:,it1)
        do n=1,nmegqwan              
          n1=bmegqwan(1,n)
          n2=bmegqwan(2,n)
          ias1=iwann(1,n1)
          ias2=iwann(1,n2)
          pos1(:)=atposc(:,ias2ia(ias1),ias2is(ias1))
          pos2(:)=atposc(:,ias2ia(ias2),ias2is(ias2))
          dvec(:)=pos2(:)+vtrc(:)-pos1(:)
          d=sqrt(dvec(1)**2+dvec(2)**2+dvec(3)**2)
          if (sum(abs(megqwan2(n,it1,:))).gt.1d-16) then
            write(150,'("  transition ",I4," between wfs : ",2I4,"  distance: ",F12.6)')&
              n,bmegqwan(1,n),bmegqwan(2,n),d
            do ig=1,ngvecme
              write(150,'("    ig : ",I4," mewf=(",2F12.6,"), |mewf|=",F12.6)')&
                ig,megqwan2(n,it1,ig),abs(megqwan2(n,it1,ig))
            enddo
          endif
        enddo !n
      endif
    enddo
  endif
  if (wproc) then
    write(150,*)
    write(150,'("Full matrix size in local basis : ",I6)')nmegqwan*ntrmegqwan
  endif
  allocate(imegqwan2(2,nmegqwan*ntrmegqwan))
  imegqwan2=0
  nmegqwan2=0
  do it1=1,ntrmegqwan
    do n=1,nmegqwan
      if (sum(abs(megqwan2(n,it1,:))).gt.1d-16) then
        nmegqwan2=nmegqwan2+1
        imegqwan2(1,nmegqwan2)=it1
        imegqwan2(2,nmegqwan2)=n
      endif
    enddo
  enddo
  if (wproc) then
    write(150,*)
    write(150,'("Reduced matrix size in local basis : ",I6)')nmegqwan2
  endif
  allocate(vcwan(nmegqwan2,nmegqwan2))
! Coulomb matrix in local basis
  vcwan=zzero
  do i=1,nmegqwan2
    do j=1,nmegqwan2
      i1=imegqwan2(1,i)
      n1=imegqwan2(2,i)
      i2=imegqwan2(1,j)
      n2=imegqwan2(2,j)
      do ig=1,ngvecme
        vcwan(i,j)=vcwan(i,j)+dconjg(megqwan2(n1,i1,ig))*megqwan2(n2,i2,ig)*&
          vcgq(ig)**2
      enddo
    enddo
  enddo
  if (wproc) call flushifc(150)
endif !wannier_chi0_chi

if (crpa) then
  allocate(imegqwan(nwann,nwann))
  imegqwan=-1
  do i=1,nmegqwan
    n1=bmegqwan(1,i)
    n2=bmegqwan(2,i)
    imegqwan(n1,n2)=i
  enddo  
endif !crpa

igq0=lr_igq0-gvecme1+1

ie1=1

allocate(chi0(ngvecme,ngvecme))
allocate(krnl(ngvecme,ngvecme))
!if (screened_w) allocate(krnl_scr(ngvecchi,ngvecchi))
allocate(ixcft(ngvec))
if (allocated(f_response)) deallocate(f_response)
allocate(f_response(nf_response,nepts,nfxca))
f_response=zzero

!if (mpi_grid_root(dims=(/dim_k,dim_b/))) call write_lr_header(qnm)

! distribute nfxca between 2-nd dimension 
bs=mpi_grid_map(nfxca,dim_b,offs=idx0)
ifxc1=idx0+1
ifxc2=idx0+bs

if (wproc) then
  write(150,*)
  write(150,'("fisrt energy point : ",I4)')ie1
  call flushifc(150)
endif
! loop over energy points
do ie=ie1,nepts
  chi0=zzero
  if (wannier_chi0_chi) chi0wan_k=zzero
  call timer_start(1,reset=.true.)
  call timer_start(2,reset=.true.)
! sum over k-points
  do ikloc=1,nkptnrloc
    if (nmegqblhloc(1,ikloc).gt.0) then
! for each k-point : sum over interband transitions
      call sumchi0(ikloc,lr_w(ie),chi0)
    endif
! for response in Wannier basis
    if (wannier_chi0_chi) call sumchi0wan_k(ikloc,lr_w(ie),chi0wan_k(1,1,ikloc))
  enddo !ikloc
  call timer_stop(2)
  call timer_start(3,reset=.true.)
! sum over k-points and band transitions
  call mpi_grid_reduce(chi0(1,1),ngvecme*ngvecme,dims=(/dim_k,dim_b/))
  chi0=chi0/nkptnr/omega
  if (wannier_chi0_chi) call mpi_grid_reduce(chi0wan_k(1,1,1), &
    nmegqwan*nmegqwan*nkptnrloc,dims=(/dim_b/))
  if (crpa.and.mpi_grid_root(dims=(/dim_k,dim_b/))) &
    call genwu(ngvecme,ie,chi0,vcgq,qnm)
! compute ch0 matrix in Wannier basis
  if (mpi_grid_root((/dim_b/)).and.wannier_chi0_chi) then
    call genchi0wan(igq0,chi0wan_k,chi0wan,chi0_GqGq_wan_full)
    if (lafm) chi0wan(:,:,:)=chi0wan(:,:,:)*2.d0    
  endif
  call timer_stop(3)
! compute response functions
  if (mpi_grid_root(dims=(/dim_k/)).and..not.crpa) then
! loop over fxc
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
      call solve_chi(ngvecme,igq0,vcgq,lr_w(ie),chi0,krnl,krnl_scr, &
        f_response(1,ie,ifxc))
      if (wannier_chi0_chi.and.ifxc.eq.1) then
        call solve_chi_wan(igq0,vcgq,lr_w(ie),nmegqwan2,imegqwan2,megqwan2,vcwan,&
          chi0wan,f_response(1,ie,ifxc))
        f_response(f_chi0_wann_full,ie,ifxc)=chi0_GqGq_wan_full
      endif
    enddo
  endif  
enddo !ie

if (mpi_grid_root(dims=(/dim_k/))) then
  call mpi_grid_reduce(f_response(1,1,1),nf_response*nepts*nfxca,dims=(/dim_b/))
! write response functions to .dat file
  if (mpi_grid_root(dims=(/dim_b/))) then
    do ifxc=1,nfxca
      call write_chi(lr_igq0,ivq0m,ifxc)
    enddo
  endif
endif

call mpi_grid_barrier(dims=(/dim_k,dim_b/))

deallocate(chi0)
deallocate(krnl)
deallocate(ixcft)
if (wannier_chi0_chi) then
  deallocate(chi0wan)
  deallocate(chi0wan_k)
  deallocate(imegqwan2)
  deallocate(megqwan2)
  deallocate(vcwan)
endif
if (crpa) then
  deallocate(imegqwan)
endif
deallocate(megqblh2)
deallocate(vcgq)

30 continue
if (wproc) then
  write(150,*)
  write(150,'("Done.")')
  call flushifc(150)
endif

 
return
end

#endif
