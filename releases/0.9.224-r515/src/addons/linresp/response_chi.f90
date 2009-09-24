#ifdef _HDF5_
subroutine response_chi(ivq0m)
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
complex(8), allocatable :: krnl_rpa(:,:)
complex(8), allocatable :: chi0w(:,:)
!complex(8), allocatable :: chi0wf(:)
!complex(8), allocatable :: chiwf(:)
real(8) fxca
complex(8), allocatable :: epsilon_(:,:,:)
complex(8), allocatable :: chi_(:,:,:)
real(8) fourpiq0
complex(8), allocatable :: mewf2(:,:,:)
complex(8), allocatable :: mewf4(:,:,:)
complex(8), allocatable :: mewf2_t(:,:,:)
complex(8), allocatable :: mewfx(:,:,:)
complex(8), allocatable :: mtrx_v(:,:)
!complex(8), allocatable :: epswf(:)

complex(8) zt1
integer ntr1,ntr2
integer, allocatable :: itr1l(:,:)
integer, allocatable :: itr2l(:,:)
integer, allocatable :: itridx(:,:)

integer nwfme
integer, allocatable :: iwfme(:,:)
logical exist
integer ntrans
integer, allocatable :: itrans(:,:)
real(8) d1
integer itrans_m
integer nnzme 
integer, allocatable :: inzme(:,:)

integer ie,ig,i,j,ig1,ig2,ie1,ie2,idx0,bs,n1,n2,it1,it2,n3,n4,n
integer i1,i2,ifxc,ifxc1,ifxc2
integer iv(3)
character*100 fname,qnm,path
logical, external :: root_cart
complex(8), external :: zdotu

! we need Bcx and magnetization from STATE.OUT
if (lrtype.eq.1) call readstate

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

call qname(ivq0m,qnm)
qnm="./"//trim(qnm)//"/"//trim(qnm)
fname=trim(qnm)//"_chi0.hdf5"
if (root_cart((/1,1,0/))) then
  call read_integer(nepts,1,trim(fname),'/parameters','nepts')
  call read_integer(lr_igq0,1,trim(fname),'/parameters','lr_igq0')
  call read_integer(gshme1,1,trim(fname),'/parameters','gshme1')
  call read_integer(gshme2,1,trim(fname),'/parameters','gshme2')
  call read_integer(gvecme1,1,trim(fname),'/parameters','gvecme1')
  call read_integer(gvecme2,1,trim(fname),'/parameters','gvecme2')
  call read_integer(ngvecme,1,trim(fname),'/parameters','ngvecme')
  call read_real8(vq0l,3,trim(fname),'/parameters','vq0l')
  call read_real8(vq0rl,3,trim(fname),'/parameters','vq0rl')
  call read_real8(vq0c,3,trim(fname),'/parameters','vq0c')
  call read_real8(vq0rc,3,trim(fname),'/parameters','vq0rc')
  if (lwannresp) then
    call read_integer(ntr1,1,trim(fname),'/wann','ntr1')
    call read_integer(ntr2,1,trim(fname),'/wann','ntr2')
    call read_integer(nwfme,1,trim(fname),'/wann','nwfme')
  endif
endif
call i_bcast_cart(comm_cart_110,nepts,1)
call i_bcast_cart(comm_cart_110,lr_igq0,1)
call i_bcast_cart(comm_cart_110,gshme1,1)
call i_bcast_cart(comm_cart_110,gshme2,1)
call i_bcast_cart(comm_cart_110,gvecme1,1)
call i_bcast_cart(comm_cart_110,gvecme2,1)
call i_bcast_cart(comm_cart_110,ngvecme,1)
call d_bcast_cart(comm_cart_110,vq0l,3)
call d_bcast_cart(comm_cart_110,vq0rl,3)
call d_bcast_cart(comm_cart_110,vq0c,3)
call d_bcast_cart(comm_cart_110,vq0rc,3)

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

! for response in Wannier basis
if (lwannresp) then
  call i_bcast_cart(comm_cart_110,ntr1,1)
  call i_bcast_cart(comm_cart_110,ntr2,1)
  call i_bcast_cart(comm_cart_110,nwfme,1)
  allocate(itr1l(3,ntr1))
  allocate(itr2l(3,ntr2))
  allocate(itridx(ntr1,ntr1))
  allocate(iwfme(2,nwfme))
  if (root_cart((/1,1,0/))) then
    call read_integer_array(itr1l,2,(/3,ntr1/),trim(fname),'/wann','itr1')
    call read_integer_array(itr2l,2,(/3,ntr2/),trim(fname),'/wann','itr2')
    call read_integer_array(itridx,2,(/ntr1,ntr1/),trim(fname),'/wann','itridx')
    call read_integer_array(iwfme,2,(/2,nwfme/),trim(fname),'/wann','iwfme')
  endif
  call i_bcast_cart(comm_cart_110,itr1l,3*ntr1)
  call i_bcast_cart(comm_cart_110,itr2l,3*ntr2)
  call i_bcast_cart(comm_cart_110,itridx,ntr1*ntr1)
  call i_bcast_cart(comm_cart_110,iwfme,2*nwfme)
  allocate(mewf2(nwfme,ntr1,ngvecme))
  allocate(mewf4(nwfme,nwfme,ntr2))
  if (root_cart((/1,1,0/))) then
    call read_real8_array(mewf2,4,(/2,nwfme,ntr1,ngvecme/), &
      trim(fname),'/wann','mewf2')
  endif
  call d_bcast_cart(comm_cart_110,mewf2,2*nwfme*ntr1*ngvecme)
! cut off some matrix elements
  inquire(file='mewf.in',exist=exist)
  if (exist) then
    open(70,file='mewf.in',form='formatted',status='old')
    read(70,*)itrans_m
    if (itrans_m.eq.1) then
      read(70,*)d1
      if (wproc) then
        write(150,'("minimal matrix element : ",F12.6)')d1
      endif
      do it1=1,ntr1
        do n=1,nwfme
          do ig=1,ngvecme
            if (abs(mewf2(n,it1,ig)).lt.d1) mewf2(n,it1,ig)=zzero
          enddo
        enddo
      enddo
    endif
    if (itrans_m.eq.2) then
      read(70,*)ntrans
      allocate(itrans(5,ntrans))
      do i=1,ntrans
        read(70,*)itrans(1:5,i)
      enddo
      allocate(mewf2_t(nwfme,ntr1,ngvecme))
      do it1=1,ntr1
        do n=1,nwfme
          n1=iwfme(1,n)
          n2=iwfme(2,n)
          do i=1,ntrans
            if (n1.eq.itrans(1,i).and.&
                n2.eq.itrans(2,i).and.&
                itr1l(1,it1).eq.itrans(3,i).and.&
                itr1l(2,it1).eq.itrans(4,i).and.&
                itr1l(3,it1).eq.itrans(5,i)) then
              mewf2_t(n,it1,:)=mewf2(n,it1,:)
            endif
            if (n2.eq.itrans(1,i).and.&
                n1.eq.itrans(2,i).and.&
                itr1l(1,it1).eq.-itrans(3,i).and.&
                itr1l(2,it1).eq.-itrans(4,i).and.&
                itr1l(3,it1).eq.-itrans(5,i)) then
              mewf2_t(n,it1,:)=mewf2(n,it1,:)
            endif
          enddo
        enddo
      enddo
      mewf2=mewf2_t
      deallocate(mewf2_t)
      deallocate(itrans)
    endif
  endif
  if (wproc) then
    write(150,*)
    write(150,'("Matrix elements in WF basis : ")')
    write(150,*)
    do it1=1,ntr1
      if (sum(abs(mewf2(:,it1,:))).gt.1d-8) then
        write(150,'("translation : ",3I4)')itr1l(:,it1)
        do n=1,nwfme
          if (sum(abs(mewf2(n,it1,:))).gt.1d-8) then
            write(150,'("  transition ",I4," between wfs : ",2I4)')&
              n,iwfme(1,n),iwfme(2,n)
            do ig=1,ngvecme
              write(150,'("    ig : ",I4," mewf=(",2F12.6,"), |mewf|=",F12.6)')&
                ig,mewf2(n,it1,ig),abs(mewf2(n,it1,ig))
            enddo
          endif
        enddo !n
      endif
    enddo
  endif
  if (wproc) then
    write(150,*)
    write(150,'("Full matrix size in local basis : ",I6)')nwfme*ntr1
  endif
  allocate(inzme(2,nwfme*ntr1))
  inzme=0
  nnzme=0
  do it1=1,ntr1
    do n=1,nwfme
      if (sum(abs(mewf2(n,it1,:))).gt.1d-8) then
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
  allocate(mtrx_v(nnzme,nnzme))
! Coulomb matrix in local basis
  mtrx_v=zzero
  do i=1,nnzme
    do j=1,nnzme
      i1=inzme(1,i)
      n1=inzme(2,i)
      i2=inzme(1,j)
      n2=inzme(2,j)
      do ig=1,ngvecchi
        vgq0c(:)=vgc(:,ig+gvecchi1-1)+vq0rc(:)
        gq0=sqrt(vgq0c(1)**2+vgq0c(2)**2+vgq0c(3)**2)
        mtrx_v(i,j)=mtrx_v(i,j)+dconjg(mewf2(n1,i1,ig))*mewf2(n2,i2,ig)*&
          fourpi/gq0**2
      enddo
    enddo
  enddo
  if (wproc) call flushifc(150)
endif !lwannresp

!if (lwannopt) then
!  allocate(mewfx(3,nwann*nwann,ntr1))
!  if (root_cart((/1,1,0/))) then
!    call read_real8_array(mewfx,4,(/2,3,nwann*nwann,ntr1/), &
!      trim(fname),'/wann','mewfx')
!  endif
!  call d_bcast_cart(comm_cart_110,mewfx,2*3*nwann*nwann*ntr1)
!!  allocate(mtrx1(nwann*nwann*ntr1,nwann*nwann*ntr1))
!!  allocate(mtrx2(nwann*nwann*ntr1,nwann*nwann*ntr1)) 
!  allocate(epswf(nepts))
!  epswf=zzero
!endif

igq0=lr_igq0-gvecchi1+1
fourpiq0=fourpi/(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)
ig1=gvecchi1-gvecme1+1
ig2=ig1+ngvecchi-1


allocate(krnl(ngvecchi,ngvecchi))
allocate(krnl_rpa(ngvecchi,ngvecchi))
allocate(ixcft(ngvec))
allocate(lr_w(nepts))
allocate(chi0w(ngvecme,ngvecme))  
allocate(chi0m(ngvecchi,ngvecchi))
allocate(chi_(7,nepts,nfxca))
allocate(epsilon_(5,nepts,nfxca))
chi_=zzero
epsilon_=zzero
lr_w=zzero

! construct RPA kernel of the matrix equation
krnl_rpa=zzero
! for charge response
if (lrtype.eq.0) then
  do ig=1,ngvecchi
! generate G+q vectors  
    vgq0c(:)=vgc(:,ig+gvecchi1-1)+vq0rc(:)
    gq0=sqrt(vgq0c(1)**2+vgq0c(2)**2+vgq0c(3)**2)
    krnl_rpa(ig,ig)=fourpi/gq0**2
  enddo !ig
endif !lrtype.eq.0
! for magnetic response
if (lrtype.eq.1) then
  call genixc(ixcft)
! contruct Ixc_{G,G'}=Ixc(G-G')
  do i=1,ngvecchi
    do j=1,ngvecchi
      iv(:)=-ivg(:,gvecchi1+i-1)+ivg(:,gvecchi1+j-1)
      krnl_rpa(i,j)=ixcft(ivgig(iv(1),iv(2),iv(3)))
    enddo
  enddo
endif !lrtype.eq.1

! distribute energy points between 1-st dimension
call idxbos(nepts,mpi_dims(1),mpi_x(1)+1,idx0,bs)
ie1=idx0+1
ie2=idx0+bs
! distribute nfxca between 2-nd dimension 
call idxbos(nfxca,mpi_dims(2),mpi_x(2)+1,idx0,bs)
ifxc1=idx0+1
ifxc2=idx0+bs

! main loop over energy points 
do ie=ie1,ie2
! read chi0(iw) from file
  if (root_cart((/0,1,0/))) then
#ifndef _PIO_
    do i=0,mpi_dims(1)-1
    do j=0,mpi_dims(3)-1
      if (mpi_x(1).eq.i.and.mpi_x(3).eq.j) then
#endif
        write(path,'("/iw/",I8.8)')ie
        call read_real8(lr_w(ie),2,trim(fname),trim(path),'w')
        call read_real8_array(chi0w,3,(/2,ngvecme,ngvecme/), &
          trim(fname),trim(path),'chi0')
        if (lwannresp) then
          call read_real8_array(chi_(5,ie,1),1,(/2/),trim(fname),trim(path),'chi0wf')
          call read_real8_array(mewf4,4,(/2,nwfme,nwfme,ntr2/), &
            trim(fname),trim(path),'mewf4')
        endif
#ifndef _PIO_      
      endif
      call barrier(comm_cart_101)
    enddo
    enddo
#endif
  endif
  call d_bcast_cart(comm_cart_010,chi0w,2*ngvecme*ngvecme)
! prepare chi0
  chi0m(1:ngvecchi,1:ngvecchi)=chi0w(ig1:ig2,ig1:ig2)
! loop over fxc
  do ifxc=ifxc1,ifxc2
    fxca=fxca0+(ifxc-1)*fxca1
! prepare fxc kernel
    krnl=krnl_rpa
    if (lrtype.eq.0) then
      do ig=1,ngvecchi
        if (fxctype.eq.1) then
          krnl(ig,ig)=krnl(ig,ig)-fxca/2.d0
        endif
        if (fxctype.eq.2) then
! generate G+q vector  
          vgq0c(:)=vgc(:,ig+gvecchi1-1)+vq0rc(:)
          gq0=sqrt(vgq0c(1)**2+vgq0c(2)**2+vgq0c(3)**2)
          krnl(ig,ig)=krnl(ig,ig)-fourpi*fxca/gq0**2
        endif
      enddo
    endif
    call solve_chi(ngvecchi,igq0,fourpiq0,chi0m,krnl,chi_(1,ie,ifxc), &
      epsilon_(1,ie,ifxc))
    if (lwannresp.and.ifxc.eq.1) then
      call solve_chi_wf(ntr1,ntr2,itridx,nwfme,nnzme,inzme,mewf2,mewf4,mtrx_v,&
        chi_(6,ie,1),chi_(7,ie,1),igq0)
    endif
  enddo !ifxc
enddo !ie
call d_reduce_cart(comm_cart_100,.false.,lr_w,2*nepts)
call d_reduce_cart(comm_cart_110,.false.,chi_,2*7*nepts*nfxca)
call d_reduce_cart(comm_cart_110,.false.,epsilon_,2*5*nepts*nfxca)
! write response functions to .dat file
if (root_cart((/1,1,0/))) then
  do ifxc=1,nfxca
    fxca=fxca0+(ifxc-1)*fxca1
    call write_chi(lr_igq0,ivq0m,chi_(1,1,ifxc),epsilon_(1,1,ifxc),fxca)
  enddo
endif

deallocate(krnl,krnl_rpa,ixcft)
deallocate(lr_w,chi0m,chi0w,chi_,epsilon_)
return
end  

#endif
