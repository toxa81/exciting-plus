#ifdef _HDF5_
subroutine genchi0(ivq0m)
use modmain
use hdf5
implicit none
! arguments
integer, intent(in) :: ivq0m(3)
! local variables
complex(8), allocatable :: chi0w(:,:)
integer i,j,ik,ie,nkptnr_,i1,i2,i3,ikloc,nspinor_
integer idx0,bs
integer igq0
complex(8) wt
character*100 path,qnm,fout,fchi0
character*8 c8
integer ierr
logical exist
integer ie1,n,n1,n2,n3,n4,jk,ist1,ist2,ig,sz1,sz2
!complex(8), allocatable :: wann_c1(:,:,:)
!complex(8), allocatable :: wann_c2(:,:,:)
complex(8), allocatable :: zv1(:),zm2(:,:,:)
integer iv(3)

integer it1,it2
real(8) vtrc(3)
complex(8) zt1,zt2,zt3
!complex(8), allocatable :: mewfx(:,:,:)
complex(8), allocatable :: pmat(:,:,:,:)

logical l1

logical lafm
logical, external :: bndint
logical, external :: wann_diel

! HDF5
integer(hid_t) h5_root_id
integer(hid_t) h5_w_id
integer(hid_t) h5_iw_id
integer(hid_t) h5_tmp_id

complex(8), external :: zdotu

call qname(ivq0m,qnm)
qnm="./"//trim(qnm)//"/"//trim(qnm)
wproc=.false.
if (mpi_grid_root((/dim_k,dim_b/))) then
  wproc=.true.
  fout=trim(qnm)//"_CHI0.OUT"
  open(150,file=trim(fout),form='formatted',status='replace')
endif

if (wproc) then
  write(150,*)
  write(150,'("Calculation of KS polarisability chi0")')
  write(150,*)
  write(150,'("Energy mesh parameters:")')
  write(150,'("  maximum energy [eV] : ", F9.4)')maxomega
  write(150,'("  energy step    [eV] : ", F9.4)')domega
  write(150,'("  eta            [eV] : ", F9.4)')lr_eta
  call flushifc(150)
endif
  
! setup energy mesh
nepts=1+maxomega/domega
allocate(lr_w(nepts))
do i=1,nepts
  lr_w(i)=dcmplx(domega*(i-1),lr_eta)/ha2ev
enddo

!read matrix elements
if (write_megq_file) then
  call readmegq(qnm)
endif


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
  ntrchi0wan=(4*lr_maxtr+1)**3
  allocate(itrchi0wan(3,ntrchi0wan))
  i=0
  do i1=-2*lr_maxtr,2*lr_maxtr
    do i2=-2*lr_maxtr,2*lr_maxtr
      do i3=-2*lr_maxtr,2*lr_maxtr
        i=i+1
        itrchi0wan(:,i)=(/i1,i2,i3/)
      enddo
    enddo
  enddo
  allocate(itridxwan(ntrmegqwan,ntrmegqwan))
  itridxwan=-1
  do n1=1,ntrmegqwan
    do n2=1,ntrmegqwan
      iv(:)=itrmegqwan(:,n1)-itrmegqwan(:,n2)
      do i=1,ntrchi0wan
        if (itrchi0wan(1,i).eq.iv(1).and.&
            itrchi0wan(2,i).eq.iv(2).and.&
            itrchi0wan(3,i).eq.iv(3)) then
          itridxwan(n1,n2)=i
        endif
      enddo
    enddo
  enddo
  allocate(chi0wan(nmegqwan,nmegqwan,ntrchi0wan))
  allocate(zv1(nmegqwan))
  allocate(zm2(nmegqwan,nmegqwan,nkptnrloc))
endif !wannier_chi0_chi
!
!if (lwannopt) then
!  allocate(pmat(3,nstsv,nstsv,nkptnr_loc))
!endif
!
igq0=lr_igq0-gvecme1+1

ie1=0
fchi0=trim(qnm)//"_chi0.hdf5"
if (mpi_grid_root((/dim_k,dim_b/))) then
  inquire(file=trim(fchi0),exist=exist)
  exist=.false.
  if (.not.exist) then
    call h5fcreate_f(trim(fchi0),H5F_ACC_TRUNC_F,h5_root_id,ierr)
    call h5gcreate_f(h5_root_id,'parameters',h5_tmp_id,ierr)
    call h5gclose_f(h5_tmp_id,ierr)
    if (wannier_chi0_chi) then
      call h5gcreate_f(h5_root_id,'wannier',h5_tmp_id,ierr)
      call h5gclose_f(h5_tmp_id,ierr)
    endif
    call h5gcreate_f(h5_root_id,'iw',h5_w_id,ierr)
    do i=1,nepts
      write(c8,'(I8.8)')i
      call h5gcreate_f(h5_w_id,c8,h5_iw_id,ierr)
      call h5gclose_f(h5_iw_id,ierr)
    enddo
    call h5gclose_f(h5_w_id,ierr)
    call h5fclose_f(h5_root_id,ierr)
    call write_integer(nepts,1,trim(fchi0),'/parameters','nepts')
    call write_integer(lr_igq0,1,trim(fchi0),'/parameters','lr_igq0')
    call write_integer(gshme1,1,trim(fchi0),'/parameters','gshme1')
    call write_integer(gshme2,1,trim(fchi0),'/parameters','gshme2')
    call write_integer(gvecme1,1,trim(fchi0),'/parameters','gvecme1')
    call write_integer(gvecme2,1,trim(fchi0),'/parameters','gvecme2')
    call write_integer(ngvecme,1,trim(fchi0),'/parameters','ngvecme')
    call write_real8(vq0l,3,trim(fchi0),'/parameters','vq0l')
    call write_real8(vq0rl,3,trim(fchi0),'/parameters','vq0rl')
    call write_real8(vq0c,3,trim(fchi0),'/parameters','vq0c')
    call write_real8(vq0rc,3,trim(fchi0),'/parameters','vq0rc')
    call write_integer(ie1,1,trim(fchi0),'/parameters','ie1')
  else
    call read_integer(ie1,1,trim(fchi0),'/parameters','ie1')
  endif
endif
call mpi_grid_bcast(ie1,dims=(/dim_k,dim_b/))
ie1=ie1+1

allocate(chi0w(ngvecme,ngvecme))

if (wproc) then
  write(150,*)
  write(150,'("Starting chi0 summation")')
  write(150,'("  first energy point : ",I4)')ie1
  call flushifc(150)
endif
! loop over energy points
do ie=ie1,nepts
  chi0w=zzero
  if (wannier_chi0_chi) zm2=zzero
  sz1=0
  sz2=0
  call timer_start(1,reset=.true.)
  call timer_start(2,reset=.true.)
! sum over k-points
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    jk=idxkq(1,ik)
    if (nmegqblh(ikloc).gt.0) then
      bs=mpi_grid_map(nmegqblh(ikloc),dim_b,offs=idx0)
      i1=idx0+1
      i2=idx0+bs
      sz1=sz1+bs
! for each k-point : sum over interband transitions
      call sumchi0(ikloc,ik,jk,i1,i2,lr_w(ie),chi0w)
    endif
! for response in Wannier basis
    if (wannier_chi0_chi) then
      sz2=sz2+nmegqblh(ikloc)*nmegqwan**2
      do i=1,nmegqblh(ikloc)
        ist1=bmegqblh(1,i,ikloc)
        ist2=bmegqblh(2,i,ikloc)
        do n=1,nmegqwan
          n1=bmegqwan(1,n)
          n2=bmegqwan(2,n)
          zv1(n)=wann_c(n1,ist1,ikloc)*dconjg(wann_c(n2,ist2,ikloc+nkptnrloc))
        enddo
        wt=(lr_occsvnr(ist1,ik)-lr_occsvnr(ist2,jk))/(lr_evalsvnr(ist1,ik) - &
            lr_evalsvnr(ist2,jk)+lr_w(ie))
        call zgerc(nmegqwan,nmegqwan,wt,zv1,1,zv1,1,zm2(1,1,ikloc),nmegqwan)
      enddo !i
    endif !wannier_chi0_chi
  enddo !ikloc
  call timer_stop(2)
  call timer_start(3,reset=.true.)
! sum over k-points and band transitions
!  if (mpi_grid_size(dim_k).gt.1.or.mpi_grid_size(dim_b).gt.1) then
  call mpi_grid_reduce(chi0w(1,1),ngvecme*ngvecme,dims=(/dim_k,dim_b/))
!  endif
  chi0w=chi0w/nkptnr/omega
! compute ch0 matrix in Wannier basis
! todo: put in a separate call 
  if (mpi_grid_root((/dim2/)).and.wannier_chi0_chi) then
! zt3 is chi0_GqGq calculated using Wannier functions expansion    
! todo: utilize second dimension of mpi grid
    zt3=zzero
! we will split matrix elements along 1-st dimention      
    bs=mpi_grid_map(nmegqwan,dim1,offs=idx0)
    n3=idx0+1
    n4=idx0+bs    
! loop over translations
    do it2=1,ntrchi0wan
      sz2=sz2+nkptnrloc*nmegqwan**2
      !zm1=zzero
      chi0wan(:,:,it2)=zzero
      do ikloc=1,nkptnrloc
        ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
! translation vector
        vtrc(:)=avec(:,1)*itrchi0wan(1,it2)+&
                avec(:,2)*itrchi0wan(2,it2)+&
                avec(:,3)*itrchi0wan(3,it2)
! phase e^{i(k+q)T}
        zt1=exp(dcmplx(0.d0,dot_product(vkcnr(:,ik)+vq0rc(:),vtrc(:))))
! chi0wan=chi0wan+e^{i(k+q)T}*zm2(k)
        call zaxpy(nmegqwan*nmegqwan,zt1,zm2(1,1,ikloc),1,chi0wan(1,1,it2),1)
      enddo !ikloc
! sum chi0wan over all k-points
      call mpi_grid_reduce(chi0wan(1,1,it2),nmegqwan*nmegqwan, &
        dims=(/dim_k/),all=.true.)
      chi0wan(:,:,it2)=chi0wan(:,:,it2)/nkptnr/omega
      if (lafm) chi0wan(:,:,it2)=chi0wan(:,:,it2)*2.d0    
! compute chi0_GqGq using the Wannier functions expansion
!  chi0=\sum_{T,T'}\sum_{n,m,n',m'} A^{*}_{nmT}(q,Gq)*chi0wan(n,m,n',m',T-T')*A_{n'm'T'}(q,Gq)
!  where A_{nmT}(q,G)=<n,0|e^{-i(G+q)x}|m,T>
      do i=1,ntrmegqwan
        do j=1,ntrmegqwan
          if (itridxwan(i,j).eq.it2) then
            sz2=sz2+ntrmegqwan*ntrmegqwan
! loop over fraction of rows
            do n2=n3,n4
! perform column times vector multiplication
!  zt1_{n,m,T,T'}=\sum_{n',m'}chi0wan(n,m,n',m',T-T')*A_{n'm'T'}(q,Gq)
              zt1=zdotu(nmegqwan,chi0wan(1,n2,it2),1,megqwan(1,i,igq0),1)
! perform vector times vector multiplication
!  zt3=zt3+\sum_{T,T'}\sum_{n,m}A^{*}_{nmT}(q,Gq)*zt1_{n,m,T,T'}
              zt3=zt3+zt1*dconjg(megqwan(n2,j,igq0))
            enddo !n2
          endif
        enddo !j
      enddo !i
    enddo !it2
! sum all rows
    call mpi_grid_reduce(zt3,dims=(/dim1/))
  endif !(root_cart((/0,1,0/)).and.wannier_chi0_chi)    
  call timer_stop(3)
! write to file
  call timer_start(4,reset=.true.)
  if (mpi_grid_root((/dim_k,dim_b/))) then
    write(path,'("/iw/",I8.8)')ie
    call write_real8(lr_w(ie),2,trim(fchi0),trim(path),'w')
    call write_real8_array(chi0w,3,(/2,ngvecme,ngvecme/), &
      trim(fchi0),trim(path),'chi0')
    if (wannier_chi0_chi) then
      call write_real8_array(zt3,1,(/2/),trim(fchi0),trim(path),'chi0wf')
      call write_real8_array(chi0wan,4,(/2,nmegqwan,nmegqwan,ntrchi0wan/), &
        trim(fchi0),trim(path),'chi0wan')
    endif
    call rewrite_integer(ie,1,trim(fchi0),'/parameters','ie1')
  endif
  call timer_stop(4)
  call timer_stop(1)
  if (wproc) then
    write(150,'("energy point ",I4," was done in ",F8.2," seconds")')ie,timer_get_value(1)
    write(150,'("  zgerc time         : ",F8.2," seconds")')timer_get_value(2)
    write(150,'("  zgerc call speed   : ",F8.2," calls/sec.")')sz1/timer_get_value(2)
    write(150,'("  zgerc memory speed : ",F8.2," Mb/sec.")')16.d0*sz1*ngvecme*ngvecme/1048576.d0/timer_get_value(2)
    write(150,'("  sync time          : ",F8.2," seconds")')timer_get_value(3)
    write(150,'("  write time         : ",F8.2," seconds")')timer_get_value(4)
    call flushifc(150)
  endif
  call mpi_grid_barrier(dims=(/dim_k,dim_b/))
enddo !ie

if (wannier_chi0_chi) then
  if (mpi_grid_root((/dim_k,dim2,0/))) then
    call write_integer(ntrchi0wan,1,trim(fchi0),'/wannier','ntrchi0wan')
    call write_integer_array(itrchi0wan,2,(/3,ntrchi0wan/),trim(fchi0),'/wannier','itrchi0wan')
    call write_integer_array(itridxwan,2,(/ntrmegqwan,ntrmegqwan/),trim(fchi0),'/wannier','itridxwan')
  endif
endif
!!
!!if (lwannopt) then
!!  allocate(mewfx(3,nwann*nwann,ntr1))
!!  mewfx=zzero
!!  do it1=1,ntr1
!!    vtrc(:)=avec(:,1)*itr1l(1,it1)+avec(:,2)*itr1l(2,it1)+avec(:,3)*itr1l(3,it1)
!!    do ikloc=1,nkptnr_loc
!!      ik=ikptnrloc(mpi_x(1),1)+ikloc-1
!!      zt1=exp(dcmplx(0.d0,-dot_product(vkcnr(:,ik)+vq0rc(:),vtrc(:))))
!!      do n1=1,nwann
!!        do n2=1,nwann
!!          do ist1=1,nstsv
!!            do ist2=1,nstsv
!!              mewfx(:,(n1-1)*nwann+n2,it1)=mewfx(:,(n1-1)*nwann+n2,it1)+&
!!                dconjg(wann_c1(n1,ist1,ikloc))*wann_c2(n2,ist2,ikloc)*zt1*&
!!                pmat(:,ist1,ist2,ikloc)*zi/(evalsvnr(ist1,ik)-evalsvnr(ist2,ik)+swidth)
!!            enddo
!!          enddo
!!        enddo
!!      enddo
!!    enddo !ikloc
!!  enddo !itr1
!!  if (root_cart((/0,1,0/)).and.mpi_dims(1).gt.1) then
!!    call d_reduce_cart(comm_cart_100,.false.,mewfx,2*3*nwann*nwann*ntr1)
!!  endif
!!  if (root_cart((/1,1,0/))) then
!!    mewfx=mewfx/nkptnr
!!    call write_real8_array(mewfx,4,(/2,3,nwann*nwann,ntr1/), &
!!      trim(fname),'/wann','mewfx')
!!  endif
!!  deallocate(mewfx)
!!endif
!
call mpi_grid_barrier(dims=(/dim_k,dim_b/))

deallocate(lr_evalsvnr)
deallocate(lr_occsvnr)
deallocate(lr_w)
deallocate(idxkq)
deallocate(nmegqblh)
deallocate(bmegqblh)
deallocate(megqblh)
deallocate(chi0w)
if (wannier_megq) then
  deallocate(wann_c)
  deallocate(itrmegqwan)
  deallocate(megqwan)
  deallocate(bmegqwan)
endif
if (wannier_chi0_chi) then
  deallocate(itrchi0wan)
  deallocate(itridxwan)
  deallocate(chi0wan)
  deallocate(zv1)
  deallocate(zm2)
endif
!if (lwannopt) then
!  deallocate(pmat)
!endif
 
return
end

#endif
