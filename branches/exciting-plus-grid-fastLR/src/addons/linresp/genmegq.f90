#ifdef _HDF5_
subroutine genmegq(ivq0m,wfsvmtloc,wfsvitloc,ngknr,igkignr,pmat)
use modmain
implicit none
! arguments
! q-vector in k-mesh coordinates
integer, intent(in) :: ivq0m(3)
complex(8), intent(in) :: wfsvmtloc(lmmaxvr,nrfmax,natmtot,nspinor,nstsv,nkptnrloc)
complex(8), intent(in) :: wfsvitloc(ngkmax,nspinor,nstsv,nkptnrloc)
integer, intent(in) :: ngknr(nkptnrloc)
integer, intent(in) :: igkignr(ngkmax,nkptnrloc)
complex(8), intent(in) :: pmat(3,nstsv,nstsv,nkptnrloc)

! allocatable arrays
integer, allocatable :: igkignr2(:)
complex(8), allocatable :: wfsvmt2(:,:,:,:,:)
complex(8), allocatable :: wfsvit2(:,:,:)

integer n1,n2,i1,i2,i3,itr,n,ist1,ist2
real(8) vtrc(3)

integer i,j,ik,jk,ig,ikstep,sz,ikloc,complete
integer ngknr2
real(8) vkq0l(3)
integer ivg1(3)
complex(8), allocatable :: gntuju(:,:,:)
integer, allocatable :: igntuju(:,:,:,:)
integer, allocatable :: ngntuju(:,:)
integer ngntujumax
complex(8) zt1
integer nkstep

integer lmaxexp
integer lmmaxexp

character*100 :: qnm,fout,fme,fu
logical l1

logical exist

! external functions
real(8), external :: r3taxi
logical, external :: wann_diel

! comment:
! the subroutine computes <psi_{n,k}|e^{-i(G+q)x}|psi_{n',k+q}> 
! 
! switch write_megq_file controls the reading and writing of ME file
! when we write ME we have two choices: write to single file or write
!  to multiple files

! maximum l for exponent expansion
lmaxexp=lmaxvr+2
lmmaxexp=(lmaxexp+1)**2

call qname(ivq0m,qnm)
qnm="./"//trim(qnm)//"/"//trim(qnm)
wproc=.false.
if (mpi_grid_root((/dim_k,dim2/))) then
  wproc=.true.
  fout=trim(qnm)//"_ME.OUT"
  open(150,file=trim(fout),form='formatted',status='replace')
endif

complete=0
fme=trim(qnm)//"_me.hdf5"
if (mpi_grid_root((/dim_k,dim2/))) then
  inquire(file=trim(fme),exist=exist)
  if (exist) then
    call read_integer(complete,1,trim(fme),'/parameters','complete')
  endif
endif
call mpi_grid_bcast(complete,dims=(/dim_k,dim2/))
if (complete.eq.1) goto 30

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
  write(150,'("Calculation of matrix elements:")')
  write(150,'("  <n,k|e^{-i(G+q)x}|n'',k+q>")')
endif

! initialize G, q and G+q vectors
call init_g_q_gq(ivq0m,lmaxexp,lmmaxexp)

! find k+q and reduce them to first BZ (this is required to utilize the 
!   periodical property of Bloch-states: |k>=|k+K>, where K is any vector 
!   of the reciprocal lattice)
allocate(idxkq(2,nkptnr))
do ik=1,nkptnr
! k+q vector
  vkq0l(:)=vklnr(:,ik)+vq0rl(:)+1d-12
! K vector
  ivg1(:)=floor(vkq0l(:))
! reduced k+q vector: k'=k+q-K
  vkq0l(:)=vkq0l(:)-ivg1(:)
! search for index of reduced k+q vector 
  do jk=1,nkptnr
    if (r3taxi(vklnr(:,jk),vkq0l).lt.epslat) then
      idxkq(1,ik)=jk
      goto 10
    endif
  enddo
  write(*,*)
  write(*,'("Error(response_me): index of reduced k+q point is not found")')
  write(*,'(" index of k-point: ",I4)')ik
  write(*,'(" K-vector: ",3I4)')ivg1
  write(*,'(" reduced k+q vector: ",3G18.10)')vkq0l
  write(*,'(" check original q-vector coordinates")')
  write(*,*)
  call pstop
10 continue
  idxkq(2,ik)=ivgig(ivg1(1),ivg1(2),ivg1(3))
enddo

! hack for q=0
!if (ivq0m(1).eq.0.and.ivq0m(2).eq.0.and.ivq0m(3).eq.0) then
!  !vq0c=(/0.0112d0,0.0223d0,0.0331d0/)
!  vq0c=(/-0.05685621179E-01,  0.05685621179E-01,  0.05685621179E-01/)
!  vq0rc=vq0c
!endif

! setup n,n' stuff
call timer_start(1,reset=.true.)
if (spinpol) then
  allocate(spinor_ud(2,nstsv,nkptnr))
  spinor_ud=0
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do j=1,nstsv
      zt1=sum(abs(evecsvloc(1:nstfv,j,ikloc)))
      if (abs(zt1).gt.1d-10) spinor_ud(1,j,ik)=1
      zt1=sum(abs(evecsvloc(nstfv+1:nstsv,j,ikloc)))
      if (abs(zt1).gt.1d-10) spinor_ud(2,j,ik)=1
    enddo
  enddo
  call mpi_grid_reduce(spinor_ud(1,1,1),2*nstsv*nkptnr,dims=(/dim_k/),all=.true.)
endif
call getmeidx(.true.)
!call mpi_grid_reduce(nmegqblhmax,dims=(/dim_k/),all=.true.,op=op_max)
allocate(nmegqblh(nkptnrloc))
allocate(bmegqblh(2,nmegqblhmax,nkptnrloc))
if (wannier_megq) then
  allocate(nmegqblhwan(nkptnrloc))
  allocate(imegqblhwan(nmegqblhmax,nkptnrloc))
  nmegqblhwan=0
  imegqblhwan=0
endif
call getmeidx(.false.)
call timer_stop(1)
if (wproc) then
  write(150,*)
  write(150,'("Maximum number of interband transitions: ",I5)')nmegqblhmax
  write(150,'("Done in ",F8.2," seconds")')timer_get_value(1)
  call flushifc(150)
endif

allocate(nmegqblhloc(2,nkptnrloc))
do ikloc=1,nkptnrloc
  nmegqblhloc(1,ikloc)=mpi_grid_map(nmegqblh(ikloc),dim_b,offs=i)
  nmegqblhloc(2,ikloc)=i
enddo
nmegqblhlocmax=maxval(nmegqblhloc)

if (wannier_megq) then
  allocate(bmegqwan(2,nwann*nwann))
  nmegqwan=0
  do n1=1,nwann
    do n2=1,nwann
      l1=.false.
! for integer occupancy numbers take only transitions between occupied and empty bands
      if (wann_diel().and.(abs(wann_occ(n1)-wann_occ(n2)).gt.1d-8)) l1=.true.
! for fractional occupancies or cRPA calculation take all transitions
      if (.not.wann_diel().or.crpa) l1=.true.
      if (l1) then
        nmegqwan=nmegqwan+1
        bmegqwan(1,nmegqwan)=n1
        bmegqwan(2,nmegqwan)=n2
      endif
    enddo
  enddo
  if (wproc) then
    write(150,*)
    write(150,'("Number of WF transitions : ",I4)')nmegqwan
  endif
! list of translations
  ntrmegqwan=(2*megqwan_maxtr+1)**3
  allocate(itrmegqwan(3,ntrmegqwan))
  i=0
  do i1=-megqwan_maxtr,megqwan_maxtr
    do i2=-megqwan_maxtr,megqwan_maxtr
      do i3=-megqwan_maxtr,megqwan_maxtr
        i=i+1
        itrmegqwan(:,i)=(/i1,i2,i3/)
      enddo
    enddo
  enddo
  allocate(megqwan(nmegqwan,ntrmegqwan,ngvecme))
  megqwan=zzero
endif

if (write_megq_file) call write_me_header(qnm)

call getmaxgnt(lmaxexp,ngntujumax)

call timer_start(1,reset=.true.)
allocate(ngntuju(nspecies,ngvecme))
allocate(igntuju(4,ngntujumax,nspecies,ngvecme))
allocate(gntuju(ngntujumax,nspecies,ngvecme))
call gengntuju(lmaxexp,ngntujumax,ngntuju,igntuju,gntuju)
call timer_stop(1)
sz=32.d0*ngntujumax*nspecies*ngvecme/1024/1024
if (wproc) then
  write(150,*)
  write(150,'("Maximum number of Gaunt-like coefficients : ",I8)')ngntujumax
  write(150,'("  size (MB) : ",I6)')sz
  write(150,'("Done in ",F8.2," seconds")')timer_get_value(1)
  call flushifc(150)
endif

sz=16.d0*ngvecme*nmegqblhlocmax*nkptnrloc/1024/1024
if (wproc) then
  write(150,*)
  write(150,'("Size of matrix elements (MB): ",I6)')sz
  call flushifc(150)
endif
allocate(megqblh(nmegqblhlocmax,ngvecme,nkptnrloc))
allocate(wfsvmt2(lmmaxvr,nrfmax,natmtot,nspinor,nstsv))
allocate(wfsvit2(ngkmax,nspinor,nstsv))
allocate(igkignr2(ngkmax))

i=0
nkstep=mpi_grid_map(nkptnr,dim_k,x=i)
if (wproc) then
  write(150,*)
  write(150,'("Starting k-point loop")')
  call flushifc(150)
endif
do ikstep=1,nkstep
  if (wproc) then
    write(150,'("k-step ",I4," out of ",I4)')ikstep,nkstep
    call flushifc(150)
  endif
! transmit wave-functions
  call timer_start(1,reset=.true.)
  call getwfkq(ikstep,wfsvmtloc,wfsvitloc,ngknr,igkignr,wfsvmt2, &
    wfsvit2,ngknr2,igkignr2)
  if (wproc) then
    write(150,'("  wave-functions are distributed")')
    call flushifc(150)
  endif
  call timer_stop(1)
! compute matrix elements  
  call timer_start(2,reset=.true.)
  if (ikstep.le.nkptnrloc) then
!    ik=mpi_grid_map(nkptnr,dim_k,loc=ikstep)
!    jk=idxkq(1,ik)
    megqblh(:,:,ikstep)=zzero
! calculate muffin-tin contribution for all combinations of n,n'
!    call timer_start(4,reset=.true.)
!    call megqblhmt(ikstep,wfsvmtloc(1,1,1,1,1,ikstep),wfsvmt2,ngntujumax,&
!      ngntuju,igntuju,gntuju)
!    call timer_stop(4)
! calculate interstitial contribution for all combinations of n,n'
!    call timer_start(5,reset=.true.)
!    call megqblhit(ikstep,ngknr(ikstep),ngknr2,igkignr(1,ikstep),igkignr2,&
!      wfsvitloc(1,1,1,ikstep),wfsvit2)
!    call timer_stop(5)
    call timer_reset(3)
    call timer_reset(4)
    call timer_reset(5)
    call genmegqblh(ikstep,ngntujumax,ngntuju,igntuju,gntuju,ngknr(ikstep),ngknr2,&
      igkignr(1,ikstep),igkignr2,wfsvmtloc(1,1,1,1,1,ikstep),wfsvmt2, &
      wfsvitloc(1,1,1,ikstep),wfsvit2)
! hack for q=0
!    if (ivq0m(1).eq.0.and.ivq0m(2).eq.0.and.ivq0m(3).eq.0) then
!      megqblh(1,:,ik1)=zzero
!      do i=1,nmegqblh(ikstep)
!        ist1=bmegqblh(1,i,ikstep)
!        ist2=bmegqblh(2,i,ikstep)
!        if (ist1.eq.ist2) megqblh(1,i,ik1)=zone
!        megqblh(1,i,ik1)=megqblh(1,i,ik1)-&
!          dot_product(vq0rc(:),pmat(:,ist1,ist2,ikstep))/&
!          (evalsvnr(ist1,ik)-evalsvnr(ist2,ik)+swidth)          
!      enddo
!    endif
! add contribution from k-point to the matrix elements of e^{-i(G+q)x} in 
!  the basis of Wannier functions
    if (wannier_megq) then
! todo: use dim2 of mpi grid for something (tranlations or G-vectors) 
      do itr=1,ntrmegqwan
        vtrc(:)=avec(:,1)*itrmegqwan(1,itr)+&
                avec(:,2)*itrmegqwan(2,itr)+&
                avec(:,3)*itrmegqwan(3,itr)
        zt1=exp(dcmplx(0.d0,-dot_product(vkcnr(:,ik)+vq0rc(:),vtrc(:))))
        do ig=1,ngvecme
          do n=1,nmegqwan
            n1=bmegqwan(1,n)
            n2=bmegqwan(2,n)
            do i1=1,nmegqblhwan(ikstep)
              i=imegqblhwan(i1,ikstep)
              ist1=bmegqblh(1,i,ikstep)
              ist2=bmegqblh(2,i,ikstep)
              megqwan(n,itr,ig)=megqwan(n,itr,ig)+dconjg(wann_c(n1,ist1,ikstep))*&
                  wann_c(n2,ist2,ikstep+nkptnrloc)*megqblh(ig,i,ikstep)*zt1
            enddo !i1
          enddo !n
        enddo !ig      
      enddo !itr
    endif !wannier
  endif !ikstep.le.nkptnrloc
  call timer_stop(2)
  if (wproc) then
    write(150,'("  time (seconds)")')
    write(150,'("    send/recv      : ",F8.2)')timer_get_value(1)
    write(150,'("    matrix elements")')
    write(150,'("      precompute MT : ",F8.2)')timer_get_value(3)
    write(150,'("      precompute IT : ",F8.2)')timer_get_value(4)
    write(150,'("      compute       : ",F8.2)')timer_get_value(5)
    write(150,'("    total for kpt  : ",F8.2)')timer_get_value(2)
    write(150,'("  speed (me/sec)   : ",F10.2)')&
      ngvecme*nmegqblh(ikstep)/timer_get_value(2)
    call flushifc(150)
  endif
enddo !ikstep

!do ikloc=1,nkptnrloc
!  write(*,*)'ikloc=',ikloc
!  do i=1,nmegqblhloc(1,ikloc)
!    write(*,*)'  ib=',i
!    write(*,*)'    me=',megqblh(i,:,ikloc)                                                                  
!  enddo                                                                                                     
!enddo  

if (wannier_megq) then
! sum over all k-points to get <n,T=0|e^{-i(G+q)x|n',T'>
  !if (mpi_grid_root((/dim2/))) then
    call mpi_grid_reduce(megqwan(1,1,1),nmegqwan*ntrmegqwan*ngvecme,&
      dims=(/dim_k,dim2/))
  !endif
  megqwan=megqwan/nkptnr
endif

if (write_megq_file) then
  if (wproc) then
    write(150,*)
    write(150,'("Writing matrix elements")')
  endif
  call timer_start(3,reset=.true.)
  call write_me(qnm,pmat)
  call timer_stop(3)
  if (wproc) write(150,'(" Done in : ",F8.2)')timer_get_value(3)
endif  

! deallocate arrays if we saved the ME file
if (write_megq_file) then
  deallocate(megqblh)
  deallocate(nmegqblh)
  deallocate(bmegqblh)
  deallocate(idxkq)
  if (wannier_megq) then
    deallocate(bmegqwan)
    deallocate(itrmegqwan)
    deallocate(megqwan)
  endif
endif

deallocate(wfsvmt2)
deallocate(wfsvit2)
deallocate(igkignr2)
deallocate(ngntuju)
deallocate(gntuju)
deallocate(igntuju)
if (spinpol) then
  deallocate(spinor_ud)
endif

if (wannier_megq) then
  deallocate(nmegqblhwan)
  deallocate(imegqblhwan)
endif


if (mpi_grid_root((/dim_k,dim2/)).and.write_megq_file) then
  complete=1
  call rewrite_integer(complete,1,trim(fme),'/parameters','complete')
endif

call mpi_grid_barrier((/dim_k,dim2/))

30 continue
if (wproc) then
  write(150,*)
  write(150,'("Done.")')
  call flushifc(150)
endif

return
end
#endif
