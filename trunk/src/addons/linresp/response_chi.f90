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
complex(8), allocatable :: chi0wf(:)
complex(8), allocatable :: lmbd(:,:)
real(8) fxca
complex(8), allocatable :: epsilon_(:,:)
complex(8), allocatable :: chi_(:,:)
real(8) fourpiq0
complex(8), allocatable :: mewf2(:,:,:)
complex(8), allocatable :: mewf4(:,:,:)
complex(8), allocatable :: mewf2_t(:,:,:)
complex(8), allocatable :: mewfx(:,:,:)
complex(8), allocatable :: mtrx1(:,:)
complex(8), allocatable :: mtrx2(:,:)
complex(8), allocatable :: epswf(:)

complex(8) zt1
integer ntr1,ntr2
integer, allocatable :: itr1l(:,:)
integer, allocatable :: itr2l(:,:)
integer, allocatable :: itridx(:,:)

integer nwf1
integer, allocatable :: iwf1(:)
integer nwf2
integer, allocatable :: iwf2(:)
integer ntr
integer, allocatable :: itr(:,:)


integer ie,ig,i,j,ig1,ig2,ie1,ie2,idx0,bs,n1,n2,it1,it2,n3,n4
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
  if (lwannresp) then
    call read_integer(ntr1,1,trim(fname),'/wann','ntr1')
    call read_integer(ntr2,1,trim(fname),'/wann','ntr2')
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
if (lwannresp) then
  call i_bcast_cart(comm_cart_110,ntr1,1)
  call i_bcast_cart(comm_cart_110,ntr2,1)
  allocate(itr1l(3,ntr1))
  allocate(itr2l(3,ntr2))
  allocate(itridx(ntr1,ntr1))
  if (root_cart((/1,1,0/))) then
    call read_integer_array(itr1l,2,(/3,ntr1/),trim(fname),'/wann','itr1')
    call read_integer_array(itr2l,2,(/3,ntr2/),trim(fname),'/wann','itr2')
  endif
  call i_bcast_cart(comm_cart_110,itr1l,3*ntr1)
  call i_bcast_cart(comm_cart_110,itr2l,3*ntr2)
  do n1=1,ntr1
    do n2=1,ntr1
      iv(:)=itr1l(:,n1)-itr1l(:,n2)
      do i=1,ntr2
        if (itr2l(1,i).eq.iv(1).and.itr2l(2,i).eq.iv(2).and.itr2l(3,i).eq.iv(3)) then
          itridx(n1,n2)=i
        endif
      enddo
    enddo
  enddo
  allocate(mewf2(nwann*nwann,ntr1,ngvecme))
  allocate(mewf4(nwann*nwann,nwann*nwann,ntr2))
  allocate(chi0wf(nepts))
  chi0wf=zzero
  if (root_cart((/1,1,0/))) then
    call read_real8_array(mewf2,4,(/2,nwann*nwann,ntr1,ngvecme/), &
      trim(fname),'/wann','mewf2')
  endif
  call d_bcast_cart(comm_cart_110,mewf2,2*nwann*nwann*ntr1*ngvecme)
! switch off some channels
  if (.false.) then
    nwf1=1
    allocate(iwf1(nwf1))
    iwf1(1)=2

    nwf2=1
    allocate(iwf2(nwf2))
    iwf2(1)=4

    ntr=2
    allocate(itr(3,ntr))
    itr(:,1)=(/-1,0,0/)
    itr(:,2)=(/0,-1,-1/)

    allocate(mewf2_t(nwann*nwann,ntr1,ngvecme))
    mewf2_t=zzero

    do it1=1,ntr1
      do i=1,ntr
        if (itr1l(1,it1).eq.itr(1,i).and.itr1l(2,it1).eq.itr(2,i).and.&
          itr1l(3,it1).eq.itr(3,i)) then
          do n1=1,nwann
            do n2=1,nwann
              do n3=1,nwf1
                do n4=1,nwf2
                  if (iwf1(n3).eq.n1.and.iwf2(n4).eq.n2) then
                    mewf2_t((n1-1)*nwann+n2,it1,:)=mewf2((n1-1)*nwann+n2,it1,:)
                  endif
                enddo
              enddo
            enddo
          enddo
        endif
      enddo
    enddo
    mewf2=mewf2_t
  endif
endif !lwannresp
if (lwannopt) then
  allocate(mewfx(3,nwann*nwann,ntr1))
  if (root_cart((/1,1,0/))) then
    call read_real8_array(mewfx,4,(/2,3,nwann*nwann,ntr1/), &
      trim(fname),'/wann','mewfx')
  endif
  call d_bcast_cart(comm_cart_110,mewfx,2*3*nwann*nwann*ntr1)
!  allocate(mtrx1(nwann*nwann*ntr1,nwann*nwann*ntr1))
!  allocate(mtrx2(nwann*nwann*ntr1,nwann*nwann*ntr1)) 
  allocate(epswf(nepts))
  epswf=zzero
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
        if (lwannresp) then
          call read_real8_array(mewf4,4,(/2,nwann*nwann,nwann*nwann,ntr2/), &
            trim(fname),trim(path),'mewf4')
        endif
#ifndef _PIO_      
      endif
      call barrier(comm_cart_101)
    enddo
    enddo
#endif
  endif
  call d_bcast_cart(comm_cart_010,lr_w(ie),2)
  call d_bcast_cart(comm_cart_010,chi0w,2*ngvecme*ngvecme)
  if (lwannresp) then
    call d_bcast_cart(comm_cart_010,mewf4,2*nwann*nwann*nwann*nwann*ntr2)
  endif  
! prepare chi0
  ig1=gvecchi1-gvecme1+1
  ig2=ig1+ngvecchi-1
  chi0m(1:ngvecchi,1:ngvecchi)=chi0w(ig1:ig2,ig1:ig2)
! solve matrix equation  
  call solve_chi(ngvecchi,igq0,fourpiq0,chi0m,krnl,chi_(1,ie), &
    epsilon_(1,ie),lmbd(1,ie))
! in Wannier basis
  if (lwannresp) then
    do it1=1,ntr1
      do it2=1,ntr1
        do n1=1,nwann*nwann
          do n2=1,nwann*nwann
            chi0wf(ie)=chi0wf(ie)+mewf4(n1,n2,itridx(it1,it2))*mewf2(n1,it1,1)*dconjg(mewf2(n2,it2,1))
          enddo
        enddo
      enddo
    enddo
  endif
  if (lwannopt) then
    do it1=1,ntr1
      do it2=1,ntr1
        do n1=1,nwann*nwann
          do n2=1,nwann*nwann
            epswf(ie)=epswf(ie)+mewf4(n1,n2,itridx(it1,it2))*mewfx(1,n1,it1)*dconjg(mewfx(1,n2,it2))
          enddo
        enddo
      enddo
    enddo 
  endif
enddo !ie

call d_reduce_cart(comm_cart_100,.false.,lr_w,2*nepts)
call d_reduce_cart(comm_cart_100,.false.,lmbd,2*ngvecchi*nepts)
call d_reduce_cart(comm_cart_100,.false.,chi_,2*4*nepts)
call d_reduce_cart(comm_cart_100,.false.,epsilon_,2*4*nepts)
if (lwannresp) then
  call d_reduce_cart(comm_cart_100,.false.,chi0wf,2*nepts)
endif
if (lwannopt) then
  call d_reduce_cart(comm_cart_100,.false.,epswf,2*nepts)
endif

if (root_cart((/1,0,0/))) then
  call write_chi(lr_igq0,ivq0m,chi_,epsilon_,lmbd)
  if (lwannresp) then
    open(130,file='chi0wann.dat',form='formatted',status='replace')
    do ie=1,nepts
      write(130,'(3G18.10)')dreal(lr_w(ie))*ha2ev,-dreal(chi0wf(ie))/ha2ev/(au2ang)**3, &
        -dimag(chi0wf(ie))/ha2ev/(au2ang)**3
    enddo
    close(130)
!    open(130,file='chi0int.dat',form='formatted',status='replace')
!    do it2=1,ntr2
!      do n1=1,nwann*nwann
!        do ie=1,nepts
!          write(130,'(2G18.10)')dreal(lr_w(ie))*ha2ev,abs(mewf4_t(n1,n1,it2,ie))
!        enddo
!        write(130,*)
!      enddo
!    enddo
!    close(130)
    do it1=1,ntr1
      write(*,*)'itr=',it1
      do n1=1,nwann
        write(*,'(10(2F12.6))')(mewf2((n1-1)*nwann+n2,it1,1),n2=1,nwann)
      enddo
    enddo
  endif
  if (lwannopt) then
    open(130,file='epswf.dat',form='formatted',status='replace')
    do ie=1,nepts
      write(130,'(3G18.10)')dreal(lr_w(ie))*ha2ev,-dreal(epswf(ie)),-dimag(epswf(ie))
    enddo
    close(130)
  endif
endif

deallocate(krnl,ixcft)
deallocate(lr_w,chi0m,chi0w,chi_,epsilon_,lmbd)
return
end  

#endif