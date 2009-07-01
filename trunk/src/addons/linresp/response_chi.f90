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
complex(8), allocatable :: chi0w(:,:)
complex(8), allocatable :: lmbd(:,:)
real(8) fxca
complex(8), allocatable :: epsilon_(:,:)
complex(8), allocatable :: chi_(:,:)
real(8) fourpiq0

integer ie,ig,i,j,ig1,ig2,ie1,ie2,idx0,bs
integer iv(3)
character*100 fname,qnm,path
logical, external :: root_cart

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

allocate(krnl(ngvecchi,ngvecchi))
allocate(ixcft(ngvec))
! construct kernel of the matrix equation
krnl=dcmplx(0.d0,0.d0)
! for charge response
if (lrtype.eq.0) then
  fxca=fxca0+fxca1*mpi_x(2)
  if (wproc) then
    write(150,*)
    write(150,'("fxc A : ",F8.4)')fxca
    write(150,*)
    write(150,'("Coulomb potential matrix elements:")')
    write(150,'("   ig        |G+q|      |v+fxc|   ")')
    write(150,'(" -------------------------------- ")')
  endif
  do ig=1,ngvecchi
! generate G+q vectors  
    vgq0c(:)=vgc(:,ig+gvecchi1-1)+vq0rc(:)
    gq0=sqrt(vgq0c(1)**2+vgq0c(2)**2+vgq0c(3)**2)
    krnl(ig,ig)=fourpi/gq0**2
    if (fxctype.eq.1) then
      krnl(ig,ig)=krnl(ig,ig)-fxca/2.d0
    endif
    if (fxctype.eq.2) then
      krnl(ig,ig)=krnl(ig,ig)-fourpi*fxca/gq0**2
    endif
    if (wproc) then
      write(150,'(1X,I4,2X,2F12.6)')ig,gq0,abs(krnl(ig,ig))
    endif
  enddo
endif
! for magnetic response
if (lrtype.eq.1) then
  if (wproc) then
    write(150,*)
    write(150,'("Ixc matrix elements:")')
    write(150,'("   ig      |Ixc(G)|   ")')
    write(150,'(" -------------------- ")')
  endif
  call genixc(ixcft)
! contruct Ixc_{G,G'}=Ixc(G-G')
  do i=1,ngvecchi
    do j=1,ngvecchi
      ig1=gvecchi1+i-1
      ig2=gvecchi1+j-1
      iv(:)=-ivg(:,ig1)+ivg(:,ig2)
      krnl(i,j)=ixcft(ivgig(iv(1),iv(2),iv(3)))
    enddo
  enddo
  if (wproc) then
    do ig=1,2*ngvecchi
      write(150,'(1X,I4,2X,F12.6)')ig,abs(ixcft(ig))
    enddo
  endif
endif
if (wproc) then
  call flushifc(150)
endif

igq0=lr_igq0-gvecchi1+1
fourpiq0=fourpi/(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)

call idxbos(nepts,mpi_dims(1),mpi_x(1)+1,idx0,bs)
ie1=idx0+1
ie2=idx0+bs

allocate(lr_w(nepts))
allocate(chi0w(ngvecme,ngvecme))  
allocate(chi0m(ngvecchi,ngvecchi))
allocate(lmbd(ngvecchi,nepts))
allocate(chi_(4,nepts))
allocate(epsilon_(4,nepts))
lr_w=zzero
lmbd=zzero
chi_=zzero
epsilon_=zzero

do ie=ie1,ie2
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
#ifndef _PIO_      
      endif
      call barrier(comm_cart_101)
    enddo
    enddo
#endif
  endif
  call d_bcast_cart(comm_cart_010,lr_w(ie),2)
  call d_bcast_cart(comm_cart_010,chi0w,2*ngvecme*ngvecme)
! prepare chi0
  ig1=gvecchi1-gvecme1+1
  ig2=ig1+ngvecchi-1
  chi0m(1:ngvecchi,1:ngvecchi)=chi0w(ig1:ig2,ig1:ig2)
! solve matrix equation  
  call solve_chi(ngvecchi,igq0,fourpiq0,chi0m,krnl,chi_(1,ie), &
    epsilon_(1,ie),lmbd(1,ie))
enddo

call d_reduce_cart(comm_cart_100,.false.,lr_w,2*nepts)
call d_reduce_cart(comm_cart_100,.false.,lmbd,2*ngvecchi*nepts)
call d_reduce_cart(comm_cart_100,.false.,chi_,2*4*nepts)
call d_reduce_cart(comm_cart_100,.false.,epsilon_,2*4*nepts)

if (root_cart((/1,0,0/))) then
  call write_chi(lr_igq0,ivq0m,chi_,epsilon_,lmbd)
endif

deallocate(krnl,ixcft)
deallocate(lr_w,chi0m,chi0w,chi_,epsilon_,lmbd)
return
end  

#endif