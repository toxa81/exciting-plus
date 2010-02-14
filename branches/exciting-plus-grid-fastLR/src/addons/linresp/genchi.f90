#ifdef _HDF5_
subroutine genchi(ivq0m)
use modmain
implicit none
integer, intent(in) :: ivq0m(3)

integer igq0
! Kohn-Sham polarisability submatrix
complex(8), allocatable :: chi0m(:,:)
! G+q vector in Cartesian coordinates
real(8) vgq0c(3)
! length of G+q vector
real(8) gq0
! Fourier-transform of Ixc(r)=Bxc(r)/m(r)
complex(8), allocatable :: ixcft(:)

! kernel of matrix equaiton
complex(8), allocatable :: krnl(:,:)
complex(8), allocatable :: krnl_scr(:,:)
complex(8), allocatable :: chi0w(:,:)
real(8) fxca
real(8) fourpiq0
complex(8), allocatable :: megqwan_t(:,:,:)
!complex(8), allocatable :: mewfx(:,:,:)
complex(8), allocatable :: vcwan(:,:)
complex(8), allocatable :: uscrn(:,:)
complex(8), allocatable :: ubare(:,:)
real(8), allocatable :: vcgq(:)

complex(8), allocatable :: chi0wan(:,:,:)

complex(8) zt1
real(8) dvec(3),pos1(3),pos2(3),vtrc(3)
integer ias1,ias2
real(8) d

logical exist
integer ntrans
integer, allocatable :: itrans(:,:)
integer itrans_m
integer nnzme 
integer, allocatable :: inzme(:,:)

integer ig,i,j,ig1,ig2,idx0,bs,n1,n2,it1,it2,n3,n4,n
integer i1,i2,ifxc,ifxc1,ifxc2
integer iv(3)
character*100 fchi,fchi0,fme,qnm,path,fw
character*3 c3
integer nwstep,iwstep,iw,nwloc
complex(8), external :: zdotu

! we need Bcx and magnetization from STATE.OUT
if (lrtype.eq.1) call readstate

call qname(ivq0m,qnm)
qnm="./"//trim(qnm)//"/"//trim(qnm)
wproc=.false.
if (mpi_grid_root((/dim_w/))) then
  wproc=.true.
  write(c3,'(I3.3)')mpi_grid_x(dim_f)
  fchi=trim(qnm)//"_CHI_"//c3//".OUT"
  open(150,file=trim(fchi),form='formatted',status='replace')
endif

if (wproc) then
  write(150,*)
  if (lrtype.eq.0) then
    write(150,'("Calculation of charge polarizability chi")')
    write(150,'("  type of fxc kernel : ")')
    if (fxctype.eq.0) write(150,'("    fxc=0 (RPA)")')
    if (fxctype.eq.1) write(150,'("    fxc=-A/2 \delta_{GG''}")')
    if (fxctype.eq.2) write(150,'("    fxc=-4*Pi*A/|G+q| \delta_{GG''}")')
  endif
  if (lrtype.eq.1) then
    write(150,'("Calculation of magnetic polarizability chi")')  
  endif
  write(150,*)
  call flushifc(150)
endif

fchi0=trim(qnm)//"_chi0.hdf5"
fme=trim(qnm)//"_me.hdf5"

! read from chi0 file
if (mpi_grid_root(dims=(/dim_w,dim_f/))) then
  call read_integer(nepts,1,trim(fchi0),'/parameters','nepts')
  call read_integer(lr_igq0,1,trim(fchi0),'/parameters','lr_igq0')
  call read_integer(gshme1,1,trim(fchi0),'/parameters','gshme1')
  call read_integer(gshme2,1,trim(fchi0),'/parameters','gshme2')
  call read_integer(gvecme1,1,trim(fchi0),'/parameters','gvecme1')
  call read_integer(gvecme2,1,trim(fchi0),'/parameters','gvecme2')
  call read_integer(ngvecme,1,trim(fchi0),'/parameters','ngvecme')
  call read_real8(vq0l,3,trim(fchi0),'/parameters','vq0l')
  call read_real8(vq0rl,3,trim(fchi0),'/parameters','vq0rl')
  call read_real8(vq0c,3,trim(fchi0),'/parameters','vq0c')
  call read_real8(vq0rc,3,trim(fchi0),'/parameters','vq0rc')
  if (wannier_chi0_chi) then
    call read_integer(ntrchi0wan,1,trim(fchi0),'/wannier','ntrchi0wan')
  endif
endif
call mpi_grid_bcast(nepts,dims=(/dim_w,dim_f/))
call mpi_grid_bcast(lr_igq0,dims=(/dim_w,dim_f/))
call mpi_grid_bcast(gshme1,dims=(/dim_w,dim_f/))
call mpi_grid_bcast(gshme2,dims=(/dim_w,dim_f/))
call mpi_grid_bcast(gvecme1,dims=(/dim_w,dim_f/))
call mpi_grid_bcast(gvecme2,dims=(/dim_w,dim_f/))
call mpi_grid_bcast(ngvecme,dims=(/dim_w,dim_f/))
call mpi_grid_bcast(vq0l(1),3,dims=(/dim_w,dim_f/))
call mpi_grid_bcast(vq0rl(1),3,dims=(/dim_w,dim_f/))
call mpi_grid_bcast(vq0c(1),3,dims=(/dim_w,dim_f/))
call mpi_grid_bcast(vq0rc(1),3,dims=(/dim_w,dim_f/))
if (wannier_chi0_chi) then
  call mpi_grid_bcast(ntrchi0wan,dims=(/dim_w,dim_f/))
endif

if (wannier_megq.and.write_megq_file) call readmegqwan(qnm)

! read arrays from chi0 file
if (wannier_chi0_chi) then
  allocate(itrchi0wan(3,ntrchi0wan))
  allocate(itridxwan(ntrmegqwan,ntrmegqwan))
  if (mpi_grid_root(dims=(/dim_w,dim_f/))) then
    call read_integer_array(itrchi0wan,2,(/3,ntrchi0wan/),trim(fchi0),'/wannier','itrchi0wan')
    call read_integer_array(itridxwan,2,(/ntrmegqwan,ntrmegqwan/),trim(fchi0),'/wannier','itridxwan')
  endif
  call mpi_grid_bcast(itrchi0wan(1,1),3*ntrchi0wan,dims=(/dim_w,dim_f/))
  call mpi_grid_bcast(itridxwan(1,1),ntrmegqwan*ntrmegqwan,dims=(/dim_w,dim_f/))
endif

if (wproc) then
  write(150,'("chi0 was calculated for ")')
  write(150,'("  G-shells  : ",I4," to ",I4)')gshme1,gshme2
  write(150,'("  G-vectors : ",I4," to ",I4)')gvecme1,gvecme2
  call flushifc(150)
endif

gshchi1=gshme1
gshchi2=gshme2
gvecchi1=gvecme1
gvecchi2=gvecme2
ngvecchi=gvecchi2-gvecchi1+1  
if (wproc) then 
  write(150,*)
  write(150,'("Minimum and maximum G-vectors for chi : ",2I4)')gvecchi1,gvecchi2
  write(150,'("Number of G-vectors : ",I4)')ngvecchi
  call flushifc(150)
endif

igq0=lr_igq0-gvecchi1+1
ig1=gvecchi1-gvecme1+1
ig2=ig1+ngvecchi-1

allocate(krnl(ngvecchi,ngvecchi))
if (screened_w) allocate(krnl_scr(ngvecchi,ngvecchi))
allocate(vcgq(ngvecchi))
allocate(ixcft(ngvec))
allocate(lr_w(nepts))
allocate(chi0w(ngvecme,ngvecme))  
allocate(chi0m(ngvecchi,ngvecchi))
allocate(f_response(nf_response,nepts,nfxca))
f_response=zzero
lr_w=zzero

! generate sqrt(4*Pi)/|G+q|  
do ig=1,ngvecchi
  vgq0c(:)=vgc(:,ig+gvecchi1-1)+vq0rc(:)
  gq0=sqrt(vgq0c(1)**2+vgq0c(2)**2+vgq0c(3)**2)
  vcgq(ig)=sqrt(fourpi)/gq0
enddo !ig

! for response in Wannier basis
if (wannier_chi0_chi) then
  allocate(chi0wan(nmegqwan,nmegqwan,ntrchi0wan))
  if (wproc) then
    write(150,'("minimal matrix element : ",F12.6)')megqwan_cutoff
  endif
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
      if (d.gt.megqwan_maxdist) megqwan(n,it1,:)=zzero
      do ig=1,ngvecme
        if (abs(megqwan(n,it1,ig)).lt.megqwan_cutoff) megqwan(n,it1,ig)=zzero
      enddo
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
      if (sum(abs(megqwan(:,it1,:))).gt.1d-8) then
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
          if (sum(abs(megqwan(n,it1,:))).gt.1d-8) then
            write(150,'("  transition ",I4," between wfs : ",2I4,"  distance: ",F12.6)')&
              n,bmegqwan(1,n),bmegqwan(2,n),d
            do ig=1,ngvecme
              write(150,'("    ig : ",I4," mewf=(",2F12.6,"), |mewf|=",F12.6)')&
                ig,megqwan(n,it1,ig),abs(megqwan(n,it1,ig))
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
  allocate(inzme(2,nmegqwan*ntrmegqwan))
  inzme=0
  nnzme=0
  do it1=1,ntrmegqwan
    do n=1,nmegqwan
      if (sum(abs(megqwan(n,it1,:))).gt.1d-8) then
        nnzme=nnzme+1
        inzme(1,nnzme)=it1
        inzme(2,nnzme)=n
      endif
    enddo
  enddo
  if (wproc) then
    write(150,*)
    write(150,'("Reduced matrix size in local basis : ",I6)')nnzme
  endif
  allocate(vcwan(nnzme,nnzme))
! Coulomb matrix in local basis
  vcwan=zzero
  do i=1,nnzme
    do j=1,nnzme
      i1=inzme(1,i)
      n1=inzme(2,i)
      i2=inzme(1,j)
      n2=inzme(2,j)
      do ig=1,ngvecchi
        vcwan(i,j)=vcwan(i,j)+dconjg(megqwan(n1,i1,ig))*megqwan(n2,i2,ig)*&
          vcgq(ig)**2
      enddo
    enddo
  enddo
  if (wproc) call flushifc(150)
endif !wannier_chi0_chi

! distribute energy points between 1-st dimension
i=0
nwstep=mpi_grid_map(nepts,dim_w,x=i)
nwloc=mpi_grid_map(nepts,dim_w)
! distribute nfxca between 2-nd dimension 
bs=mpi_grid_map(nfxca,dim_f,offs=idx0)
ifxc1=idx0+1
ifxc2=idx0+bs
! main loop over energy points 
do iwstep=1,nwstep
! read chi0(iw) from file
  if (mpi_grid_root(dims=(/dim_f/))) then
    do i=0,mpi_grid_size(dim_w)-1
    do j=0,mpi_grid_size(dim_q)-1
      if (mpi_grid_x(dim_w).eq.i.and.mpi_grid_x(dim_q).eq.j) then
        if (iwstep.le.nwloc) then
          iw=mpi_grid_map(nepts,dim_w,loc=iwstep)
          write(path,'("/iw/",I8.8)')iw
          call read_real8(lr_w(iw),2,trim(fchi0),trim(path),'w')
          call read_real8_array(chi0w,3,(/2,ngvecme,ngvecme/), &
            trim(fchi0),trim(path),'chi0')
          if (wannier_chi0_chi) then
            call read_real8_array(f_response(f_chi0_wann_full,iw,1),1,(/2/), &
              trim(fchi0),trim(path),'chi0wf')
            call read_real8_array(chi0wan,4,(/2,nmegqwan,nmegqwan,ntrchi0wan/), &
              trim(fchi0),trim(path),'chi0wan')
          endif
        endif
      endif
      if (.not.parallel_read) call mpi_grid_barrier(dims=(/dim_w,dim_q/))
    enddo !j
    enddo !i
  endif
  if (iwstep.le.nwloc) then
    iw=mpi_grid_map(nepts,dim_w,loc=iwstep)
    call mpi_grid_bcast(chi0w(1,1),ngvecme*ngvecme,dims=(/dim_f/))
! prepare chi0
    chi0m(1:ngvecchi,1:ngvecchi)=chi0w(ig1:ig2,ig1:ig2)
! loop over fxc
    do ifxc=ifxc1,ifxc2
      fxca=fxca0+(ifxc-1)*fxca1
! prepare fxc kernel
      krnl=zzero
      if (lrtype.eq.0) then
        do ig=1,ngvecchi
          if (fxctype.eq.1) then
            krnl(ig,ig)=krnl(ig,ig)-fxca/2.d0
          endif
          if (fxctype.eq.2) then
            krnl(ig,ig)=krnl(ig,ig)-fxca*vcgq(ig)**2
          endif
        enddo
      endif
      call solve_chi(igq0,vcgq,lr_w(iw),chi0m,krnl,krnl_scr,f_response(1,iw,ifxc))
      if (wannier_chi0_chi.and.ifxc.eq.1) then
        call solve_chi_wf(igq0,vcgq,lr_w(iw),nnzme,inzme,vcwan,chi0wan,f_response(1,iw,ifxc))
      endif
      if (screened_w.and.iw.eq.1.and.ifxc.eq.1) then
        if (ngvecchi.gt.10) then
          n1=10
        else
          n1=ngvecchi
        endif
        fw=trim(qnm)//"_W.txt"
        open(170,file=trim(fw),status='replace',form='formatted')
        write(170,'("Screened W matrix")')
        write(170,'("real part")')
        do i=1,n1
          write(170,'(100F12.6)')(dreal(krnl_scr(i,j)),j=1,n1)
        enddo
        write(170,'("imag part")')
        do i=1,n1
          write(170,'(100F12.6)')(dimag(krnl_scr(i,j)),j=1,n1)
        enddo
        close(170)
      endif !screened_w.and.iw.eq.1.and.ifxc.eq.1
    enddo !ifxc
  endif !iwstep.le.nwloc
enddo !ie

call mpi_grid_reduce(lr_w(1),nepts,dims=(/dim_w/))
call mpi_grid_reduce(f_response(1,1,1),nf_response*nepts*nfxca,dims=(/dim_w,dim_f/))
! write response functions to .dat file
if (mpi_grid_root(dims=(/dim_w,dim_f/))) then
  do ifxc=1,nfxca
    call write_chi(lr_igq0,ivq0m,ifxc)
  enddo
endif

deallocate(krnl,ixcft)
if (screened_w) deallocate(krnl_scr)
deallocate(lr_w,chi0m,chi0w,f_response)
deallocate(vcgq)

if (wannier_megq) then
  deallocate(itrmegqwan)
  deallocate(megqwan)
  deallocate(bmegqwan)
endif
if (wannier_chi0_chi) then
  deallocate(itrchi0wan)
  deallocate(itridxwan)
  deallocate(chi0wan)
  deallocate(vcwan)
endif
if (crpa) then
  deallocate(imegqwan)
endif
return
end  
#endif
