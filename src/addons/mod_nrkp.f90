module mod_nrkp
use mod_wannier

integer, allocatable :: ngknr(:)
integer, allocatable :: igkignr(:,:)
real(8), allocatable :: vgklnr(:,:,:)
real(8), allocatable :: vgkcnr(:,:,:)
real(8), allocatable :: gknr(:,:)
real(8), allocatable :: tpgknr(:,:,:)
complex(8), allocatable :: sfacgknr(:,:,:)
complex(8), allocatable :: ylmgknr(:,:,:)

complex(8), allocatable :: wfsvmtnrloc(:,:,:,:,:,:)
complex(8), allocatable :: wfsvitnrloc(:,:,:,:)
complex(8), allocatable :: wanncnrloc(:,:,:)
complex(8), allocatable :: pmatnrloc(:,:,:,:)

real(8), allocatable :: evalsvnr(:,:)
real(8), allocatable :: occsvnr(:,:)
integer, allocatable :: spinor_ud(:,:,:)

contains

subroutine gengknr
use modmain
implicit none
integer ik,ikloc,ig
! generate G+k vectors for entire BZ (this is required to compute 
!   wave-functions at each k-point)
if (allocated(vgklnr)) deallocate(vgklnr)
allocate(vgklnr(3,ngkmax,nkptnrloc))
if (allocated(vgkcnr)) deallocate(vgkcnr)
allocate(vgkcnr(3,ngkmax,nkptnrloc))
if (allocated(gknr)) deallocate(gknr)
allocate(gknr(ngkmax,nkptnrloc))
if (allocated(tpgknr)) deallocate(tpgknr)
allocate(tpgknr(2,ngkmax,nkptnrloc))
if (allocated(ngknr)) deallocate(ngknr)
allocate(ngknr(nkptnrloc))
if (allocated(sfacgknr)) deallocate(sfacgknr)
allocate(sfacgknr(ngkmax,natmtot,nkptnrloc))
if (allocated(igkignr)) deallocate(igkignr)
allocate(igkignr(ngkmax,nkptnrloc))
if (allocated(ylmgknr)) deallocate(ylmgknr)
allocate(ylmgknr(lmmaxvr,ngkmax,nkptnrloc))
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  call gengpvec(vklnr(1,ik),vkcnr(1,ik),ngknr(ikloc),igkignr(1,ikloc), &
    vgklnr(1,1,ikloc),vgkcnr(1,1,ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc))
  call gensfacgp(ngknr(ikloc),vgkcnr(1,1,ikloc),ngkmax,sfacgknr(1,1,ikloc))
  do ig=1,ngknr(ikloc)
    call genylm(lmaxvr,tpgknr(1,ig,ikloc),ylmgknr(1,ig,ikloc))
  enddo
enddo
return
end subroutine

subroutine genwfnr(fout,lpmat)
use modmain
implicit none
integer, intent(in) :: fout
logical, intent(in) :: lpmat
integer ik,ikloc,n,j,ik1,isym,i,ierr
complex(8), allocatable :: apwalm(:,:,:,:)
real(8) w2,t1,sz
logical, external :: wann_diel
complex(8), allocatable :: evecfvnrloc(:,:,:,:)
complex(8), allocatable :: evecsvnrloc(:,:,:)

! get energies of states in reduced part of BZ
call timer_start(3,reset=.true.)
if (wproc.and.fout.gt.0) then
  write(fout,*)
  write(fout,'("Reading energies of states")')
  if (fout.ne.6) call flushifc(fout)
endif
if (mpi_grid_root()) then
! read from IBZ
  do ik=1,nkpt
    call getevalsv(vkl(1,ik),evalsv(1,ik))
  enddo
endif
call mpi_grid_bcast(evalsv(1,1),nstsv*nkpt)
if (allocated(evalsvnr)) deallocate(evalsvnr)
allocate(evalsvnr(nstsv,nkptnr))
evalsvnr=0.d0
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  call findkpt(vklnr(1,ik),isym,ik1) 
  evalsvnr(:,ik)=evalsv(:,ik1)
enddo
call timer_stop(3)
if (wproc.and.fout.gt.0) then
  write(fout,'("Done in ",F8.2," seconds")')timer_get_value(3)
  call timestamp(fout)
  if (fout.ne.6) call flushifc(fout)
endif
!if (wproc.and.fout.gt.0) then
!  write(fout,*)
!  write(fout,'("non-reduced to reduced k-point mapping")')
!  do ik=1,nkptnr
!    call findkpt(vklnr(1,ik),isym,ik1)
!    write(fout,'(" ik : ",I4,", vklnr(ik) : ",3F10.6,", jk : ",I4,&
!      &", vkl(jk) : ",3F10.6)')ik,vklnr(:,ik),ik1,vkl(:,ik1)
!  enddo
!endif
call gengknr
if (wproc.and.fout.gt.0) then
! eigen-vectors
  sz=dble(nmatmax*nstfv*nspnfv)+dble(nstsv*nstsv)
! wave-functions
  sz=sz+dble(lmmaxvr*nufrmax*natmtot*nstsv*nspinor)+dble(ngkmax*nstsv*nspinor)
! wannier functions
  if (wannier) then
    sz=sz+dble(nwantot*nstsv)+dble(lmmaxvr*nufrmax*natmtot*nspinor*nwantot)+&
      dble(ngkmax*nspinor*nwantot)
  endif
! momentum operator matrix
  if (lpmat) then
    sz=sz+dble(3*nstsv*nstsv)
  endif
  sz=16.d0*sz*nkptnrloc/1024/1024
  write(fout,*)
  write(fout,'("Size of wave-function arrays (MB) : ",F10.2)')sz
  write(fout,*)
  write(fout,'("Reading eigen-vectors")')
  if (fout.ne.6) call flushifc(fout)
endif
call mpi_grid_barrier()
if (allocated(wfsvmtnrloc)) deallocate(wfsvmtnrloc)
allocate(wfsvmtnrloc(lmmaxvr,nufrmax,natmtot,nspinor,nstsv,nkptnrloc))
if (allocated(wfsvitnrloc)) deallocate(wfsvitnrloc)
allocate(wfsvitnrloc(ngkmax,nspinor,nstsv,nkptnrloc))
allocate(evecfvnrloc(nmatmax,nstfv,nspnfv,nkptnrloc))
allocate(evecsvnrloc(nstsv,nstsv,nkptnrloc))
if (lpmat) then
  if (allocated(pmatnrloc)) deallocate(pmatnrloc)
  allocate(pmatnrloc(3,nstsv,nstsv,nkptnrloc))
endif
if (wannier) then
  if (allocated(wanncnrloc)) deallocate(wanncnrloc)
  allocate(wanncnrloc(nwantot,nstsv,nkptnrloc))
  if (allocated(wann_unkmt)) deallocate(wann_unkmt)
  allocate(wann_unkmt(lmmaxvr,nufrmax,natmtot,nspinor,nwantot,nkptnrloc))
  if (allocated(wann_unkit)) deallocate(wann_unkit)
  allocate(wann_unkit(ngkmax,nspinor,nwantot,nkptnrloc))
endif
call timer_start(1,reset=.true.)
! read eigen-vectors
if (mpi_grid_side(dims=(/dim_k/))) then
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    call getevecfv(vklnr(1,ik),vgklnr(1,1,ikloc),evecfvnrloc(1,1,1,ikloc))
    call getevecsv(vklnr(1,ik),evecsvnrloc(1,1,ikloc))
  enddo !ikloc
endif !mpi_grid_side(dims=(/dim_k/)
call mpi_grid_barrier
call mpi_grid_bcast(evecfvnrloc(1,1,1,1),nmatmax*nstfv*nspnfv*nkptnrloc,&
  dims=ortdims((/dim_k/)))
call mpi_grid_bcast(evecsvnrloc(1,1,1),nstsv*nstsv*nkptnrloc,&
  dims=ortdims((/dim_k/)))
call timer_stop(1)
if (wproc.and.fout.gt.0) then
  write(fout,'("Done in ",F8.2," seconds")')timer_get_value(1)
  if (fout.ne.6) call flushifc(fout)
endif
call timer_start(1,reset=.true.)
if (wproc.and.fout.gt.0) then
  write(fout,*)
  write(fout,'("Generating wave-functions")')
  if (fout.ne.6) call flushifc(fout)
endif
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
! transform eigen-vectors
wfsvmtnrloc=zzero
wfsvitnrloc=zzero
if (ldisentangle) tevecsv=.true.
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
! get apw coeffs 
  call match(ngknr(ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc),        &
    sfacgknr(1,1,ikloc),apwalm)
! generate wave functions in muffin-tins
  call genwfsvmt(lmaxvr,lmmaxvr,ngknr(ikloc),evecfvnrloc(1,1,1,ikloc), &
    evecsvnrloc(1,1,ikloc),apwalm,wfsvmtnrloc(1,1,1,1,1,ikloc))
  if (wannier) then
    call genwann_c(ik,vkcnr(:,ik),evalsvnr(1,ik),wfsvmtnrloc(1,1,1,1,1,ikloc),&
      wanncnrloc(1,1,ikloc),ierr)   
    if (ldisentangle) then
! disentangle bands
      call disentangle(evalsvnr(1,ik),wanncnrloc(1,1,ikloc),&
        evecsvnrloc(1,1,ikloc))
      call genwfsvmt(lmaxvr,lmmaxvr,ngknr(ikloc),evecfvnrloc(1,1,1,ikloc), &
        evecsvnrloc(1,1,ikloc),apwalm,wfsvmtnrloc(1,1,1,1,1,ikloc))
    endif
  endif
! generate wave functions in interstitial
  call genwfsvit(ngknr(ikloc),evecfvnrloc(1,1,1,ikloc), &
    evecsvnrloc(1,1,ikloc),wfsvitnrloc(1,1,1,ikloc))
  if (lpmat) then
    call genpmat(ngknr(ikloc),igkignr(1,ikloc),vgkcnr(1,1,ikloc),&
      apwalm,evecfvnrloc(1,1,1,ikloc),evecsvnrloc(1,1,ikloc),&
      pmatnrloc(1,1,1,ikloc))
  endif    
enddo !ikloc
deallocate(apwalm)
call timer_stop(1)
if (wproc.and.fout.gt.0) then
  write(fout,'("Done in ",F8.2," seconds")')timer_get_value(1)
  call timestamp(fout)
  if (fout.ne.6) call flushifc(fout)
endif
! generate Wannier function expansion coefficients
if (wannier) then
  call timer_start(1,reset=.true.)
  if (allocated(wann_unkmt)) deallocate(wann_unkmt)
  allocate(wann_unkmt(lmmaxvr,nufrmax,natmtot,nspinor,nwantot,nkptnrloc))
  if (allocated(wann_unkit)) deallocate(wann_unkit)
  allocate(wann_unkit(ngkmax,nspinor,nwantot,nkptnrloc))
  wann_unkmt=zzero
  wann_unkit=zzero
  if (wproc.and.fout.gt.0) then
    write(fout,*)
    write(fout,'("Generating Wannier functions")')
    if (fout.ne.6) call flushifc(fout)
  endif !wproc
  do ikloc=1,nkptnrloc
    do n=1,nwantot
      do j=1,nstsv
        wann_unkmt(:,:,:,:,n,ikloc)=wann_unkmt(:,:,:,:,n,ikloc) + &
          wfsvmtnrloc(:,:,:,:,j,ikloc)*wanncnrloc(n,j,ikloc)
        wann_unkit(:,:,n,ikloc)=wann_unkit(:,:,n,ikloc) + &
          wfsvitnrloc(:,:,j,ikloc)*wanncnrloc(n,j,ikloc)
      enddo
    enddo
  enddo !ikloc
  call timer_stop(1)
endif !wannier
! after optinal band disentanglement we can finally synchronize all eigen-values
!   and compute band occupation numbers 
call mpi_grid_reduce(evalsvnr(1,1),nstsv*nkptnr,dims=(/dim_k/),all=.true.)
if (allocated(occsvnr)) deallocate(occsvnr)
allocate(occsvnr(nstsv,nkptnr))
call occupy2(nkptnr,wkptnr,evalsvnr,occsvnr)
if (mpi_grid_root()) then
  open(180,file='EIGVALNR.OUT',form='formatted',status='replace')
  write(180,'(I6," : nkptnr")') nkptnr
  write(180,'(I6," : nstsv")') nstsv
  do ik=1,nkptnr
    write(180,*)
    write(180,'(I6,4G18.10," : k-point, vkl, wkpt")') ik,vklnr(:,ik),wkptnr(ik)
    write(180,'(" (state, eigenvalue and occupancy below)")')
    do i=1,nstsv
      write(180,'(I6,2G18.10)')i,evalsvnr(i,ik),occsvnr(i,ik)
    end do
    write(180,*)
  end do
  close(180)
endif
if (wannier) then
  call timer_start(1)
! calculate Wannier function occupancies 
  wann_occ=0.d0
  wann_ene=0.d0
  do n=1,nwantot
    do ikloc=1,nkptnrloc
      ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
      do j=1,nstsv
        w2=dreal(dconjg(wanncnrloc(n,j,ikloc))*wanncnrloc(n,j,ikloc))
        wann_occ(n)=wann_occ(n)+w2*occsvnr(j,ik)/nkptnr
        wann_ene(n)=wann_ene(n)+w2*evalsvnr(j,ik)/nkptnr
      enddo
    enddo
  enddo
  call mpi_grid_reduce(wann_occ(1),nwantot,dims=(/dim_k/),all=.true.)
  call mpi_grid_reduce(wann_ene(1),nwantot,dims=(/dim_k/),all=.true.)
  if (wproc.and.fout.gt.0) then
    write(fout,'("  Wannier function occupancy and energy : ")')
    do n=1,nwantot
      write(fout,'("    n : ",I4,"  occ, ene : ",F8.6,2X,G18.10)')n,&
        wann_occ(n),wann_ene(n)
    enddo
  endif
  call timer_stop(1)
  if (wproc.and.fout.gt.0) then
    write(fout,'("  Dielectric Wannier functions : ",L1)')wann_diel()
    write(fout,*)
    write(fout,'("Done in ",F8.2," seconds")')timer_get_value(1)
    call timestamp(fout)
    if (fout.ne.6) call flushifc(fout)
  endif
endif
if (spinpol) then
  if (allocated(spinor_ud)) deallocate(spinor_ud)
  allocate(spinor_ud(2,nstsv,nkptnr))
  spinor_ud=0
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do j=1,nstsv
      t1=sum(abs(evecsvnrloc(1:nstfv,j,ikloc)))
      if (t1.gt.1d-10) spinor_ud(1,j,ik)=1
      t1=sum(abs(evecsvnrloc(nstfv+1:nstsv,j,ikloc)))
      if (t1.gt.1d-10) spinor_ud(2,j,ik)=1
    enddo
  enddo
  call mpi_grid_reduce(spinor_ud(1,1,1),2*nstsv*nkptnr,dims=(/dim_k/),all=.true.)
endif  
if (wproc.and.fout.gt.0) then
  write(fout,'("Done.")')
endif
deallocate(evecfvnrloc,evecsvnrloc)
end subroutine

end module
