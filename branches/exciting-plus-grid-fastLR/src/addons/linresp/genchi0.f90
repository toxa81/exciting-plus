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
complex(8) :: chi0_GqGq_wan_full


integer i,ik,ie,i1,i2,i3,ikloc
integer igq0
character*100 path,qnm,fout,fchi0,fu
logical exist
integer ie1,n1,n2,jk
integer iv(3)

real(8), allocatable :: vcgq(:)
real(8) vgq0c(3)
integer ig
real(8) gq0

logical lafm
logical, external :: bndint
logical, external :: wann_diel

complex(8), external :: zdotu

call qname(ivq0m,qnm)
qnm="./"//trim(qnm)//"/"//trim(qnm)
wproc=.false.
if (mpi_grid_root((/dim_k,dim_b/))) then
  wproc=.true.
  fout=trim(qnm)//"_CHI0.OUT"
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
nepts=1+int(maxomega/domega)
allocate(lr_w(nepts))
do i=1,nepts
  lr_w(i)=dcmplx(domega*(i-1),lr_eta)/ha2ev
enddo

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
  ntrchi0wan=(4*megqwan_maxtr+1)**3
  allocate(itrchi0wan(3,ntrchi0wan))
  i=0
  n1=2*megqwan_maxtr
  do i1=-n1,n1
    do i2=-n1,n1
      do i3=-n1,n1
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
  allocate(chi0wan_k(nmegqwan,nmegqwan,nkptnrloc))
endif !wannier_chi0_chi

if (crpa) then
  allocate(vcgq(ngvecme))
  do ig=1,ngvecme
! generate G+q vectors  
    vgq0c(:)=vgc(:,ig+gvecme1-1)+vq0rc(:)
    gq0=sqrt(vgq0c(1)**2+vgq0c(2)**2+vgq0c(3)**2)
    vcgq(ig)=sqrt(fourpi)/gq0
  enddo !ig
  allocate(imegqwan(nwann,nwann))
  imegqwan=-1
  do i=1,nmegqwan
    n1=bmegqwan(1,i)
    n2=bmegqwan(2,i)
    imegqwan(n1,n2)=i
  enddo  
endif

igq0=lr_igq0-gvecme1+1

ie1=0
fchi0=trim(qnm)//"_chi0.hdf5"
if (mpi_grid_root((/dim_k,dim_b/)).and.write_chi0_file) then
  inquire(file=trim(fchi0),exist=exist)
  if (.not.exist) then
    call write_chi0_header(qnm)
  else
    call read_integer(ie1,1,trim(fchi0),'/parameters','ie1')
  endif
endif
call mpi_grid_bcast(ie1,dims=(/dim_k,dim_b/))
ie1=ie1+1

allocate(chi0(ngvecme,ngvecme))

if (wproc) then
  write(150,*)
  write(150,'("Starting chi0 summation")')
  write(150,'("  first energy point : ",I4)')ie1
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
  if (crpa.and.ie.eq.1.and.mpi_grid_root(dims=(/dim_k,dim_b/))) &
    call genwu(ngvecme,chi0,vcgq,qnm)
! compute ch0 matrix in Wannier basis
  if (mpi_grid_root((/dim_b/)).and.wannier_chi0_chi) then
    call genchi0wan(igq0,chi0wan_k,chi0wan,chi0_GqGq_wan_full)
    if (lafm) chi0wan(:,:,:)=chi0wan(:,:,:)*2.d0    
  endif
  call timer_stop(3)
! write to file
  call timer_start(4,reset=.true.)
  if (mpi_grid_root((/dim_k,dim_b/)).and.write_chi0_file) then
    write(path,'("/iw/",I8.8)')ie
    call write_real8(lr_w(ie),2,trim(fchi0),trim(path),'w')
    call write_real8_array(chi0,3,(/2,ngvecme,ngvecme/), &
      trim(fchi0),trim(path),'chi0')
    if (wannier_chi0_chi) then
      call write_real8_array(chi0_GqGq_wan_full,1,(/2/),trim(fchi0),trim(path),'chi0wf')
      call write_real8_array(chi0wan,4,(/2,nmegqwan,nmegqwan,ntrchi0wan/), &
        trim(fchi0),trim(path),'chi0wan')
    endif
    call rewrite_integer(ie,1,trim(fchi0),'/parameters','ie1')
  endif
  call timer_stop(4)
  call timer_stop(1)
  !if (wproc) then
    !write(150,'("energy point ",I4," was done in ",F8.2," seconds")')ie,timer_get_value(1)
    !write(150,'("  zgerc time         : ",F8.2," seconds")')timer_get_value(2)
    !write(150,'("  zgerc call speed   : ",F8.2," calls/sec.")')sz1/timer_get_value(2)
    !write(150,'("  zgerc memory speed : ",F8.2," Mb/sec.")')16.d0*sz1*ngvecme*ngvecme/1048576.d0/timer_get_value(2)
    !write(150,'("  sync time          : ",F8.2," seconds")')timer_get_value(3)
    !write(150,'("  write time         : ",F8.2," seconds")')timer_get_value(4)
    !call flushifc(150)
  !endif
  call mpi_grid_barrier(dims=(/dim_k,dim_b/))
enddo !ie

if (wannier_chi0_chi.and.write_chi0_file) then
  if (mpi_grid_root((/dim_k,dim2/))) then
    call write_integer(ntrchi0wan,1,trim(fchi0),'/wannier','ntrchi0wan')
    call write_integer_array(itrchi0wan,2,(/3,ntrchi0wan/),trim(fchi0),'/wannier','itrchi0wan')
    call write_integer_array(itridxwan,2,(/ntrmegqwan,ntrmegqwan/),trim(fchi0),'/wannier','itridxwan')
  endif
endif

call mpi_grid_barrier(dims=(/dim_k,dim_b/))

deallocate(lr_w)
deallocate(idxkq)
deallocate(nmegqblh)
deallocate(bmegqblh)
deallocate(megqblh)
deallocate(chi0)
if (wannier_megq.and.write_megq_file) then
  deallocate(itrmegqwan)
  deallocate(megqwan)
  deallocate(bmegqwan)
endif
if (wannier_chi0_chi) then
  deallocate(itrchi0wan)
  deallocate(itridxwan)
  deallocate(chi0wan)
endif
if (crpa) then
  deallocate(vcgq)
  deallocate(imegqwan)
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
