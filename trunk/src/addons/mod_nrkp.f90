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
!
!
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
allocate(ylmgknr(lmmaxapw,ngkmax,nkptnrloc))
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  call gengpvec(vklnr(1,ik),vkcnr(1,ik),ngknr(ikloc),igkignr(1,ikloc), &
    &vgklnr(1,1,ikloc),vgkcnr(1,1,ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc))
  call gensfacgp(ngknr(ikloc),vgkcnr(1,1,ikloc),ngkmax,sfacgknr(1,1,ikloc))
  do ig=1,ngknr(ikloc)
    call genylm(lmaxapw,tpgknr(1,ig,ikloc),ylmgknr(1,ig,ikloc))
  enddo
enddo
return
end subroutine

!subroutine gen_k_sym
!use modmain
!implicit none
!integer nk,ik,isym,i,j
!real(8) s(3,3),v1(3)
!logical lfound
!integer, allocatable :: k0(:)
!integer, allocatable :: nk0(:)
!real(8), allocatable :: vkl_(:,:)
!real(8), allocatable :: wk_(:)
!
!allocate(vkl_(3,nsymcrys*nkpt))
!allocate(k0(nsymcrys*nkpt))
!allocate(nk0(nkpt))
!allocate(wk_(nsymcrys*nkpt))
!nk=0
!nk0=0
!do ik=1,nkpt
!  do isym=1,nsymcrys
!    s(:,:)=dble(symlat(:,:,lsplsymc(isym)))
!    call r3mtv(s,vkl(1,ik),v1)
!    lfound=.false.
!    do i=1,nk
!      if (sum(abs(v1(:)-vkl_(:,i))).lt.1d-8) lfound=.true.
!    enddo
!    if (.not.lfound) then
!      nk=nk+1
!      vkl_(:,nk)=v1
!      k0(nk)=ik
!      nk0(ik)=nk0(ik)+1
!    endif
!  enddo
!enddo
!wk_=0.d0
!do ik=1,nk
!  wk_(ik)=wkpt(k0(ik))/nk0(k0(ik))
!enddo
!
!deallocate(vklnr,vkcnr,wkptnr)
!nkptnr=nk
!nkptnrloc=mpi_grid_map(nkptnr,dim_k)
!allocate(vklnr(3,nkptnr),vkcnr(3,nkptnr),wkptnr(nkptnr))
!
!vklnr(:,1:nkptnr)=vkl_(:,1:nkptnr)
!wkptnr(1:nkptnr)=wk_(1:nkptnr)
!do ik=1,nkptnr
!  call r3mv(bvec,vklnr(1,ik),vkcnr(1,ik))
!enddo
!!do ik=1,nk
!!  write(*,*)vkl_(:,ik)," -- ",k0(ik),wk_(ik)
!!enddo
!deallocate(vkl_,k0,nk0,wk_)
!return
!end subroutine


subroutine wancnr_transform(umtrx)
use modmain
use mod_wannier
complex(8), intent(in) :: umtrx(nwantot,nwantot,nkptnrloc)
!
integer ikloc,ik,n,m,j
complex(8) zt1
real(8) w2
complex(8), allocatable :: wanc(:,:)
!
if (ldisentangle) then
  write(*,'("Error(wancnr_transform): disentanglement is not implemented here")')
  call pstop
endif
if (allocated(wann_unkmt)) deallocate(wann_unkmt)
allocate(wann_unkmt(lmmaxapw,nufrmax,natmtot,nspinor,nwantot,nkptnrloc))
if (allocated(wann_unkit)) deallocate(wann_unkit)
allocate(wann_unkit(ngkmax,nspinor,nwantot,nkptnrloc))
wann_unkmt=zzero
wann_unkit=zzero
allocate(wanc(nwantot,nstsv))
do ikloc=1,nkptnrloc
  do n=1,nwantot
    do m=1,nwantot
      zt1=zzero
      do j=1,nwantot
        zt1=zt1+umtrx(j,m,ikloc)*dconjg(umtrx(j,n,ikloc))
      enddo
      if (n.eq.m) zt1=zt1-zone
      if (abs(zt1).gt.1d-10) then
        write(*,'("Error(wancnr_transform): umtrx in not hermitian")')
        call pstop
      endif
    enddo
  enddo  
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc) 
  call wan_gencsv(lmmaxapw,vkcnr(1,ik),evalsvnr(1,ik),&
        &wfsvmtnrloc(1,1,1,1,1,ikloc),wanncnrloc(1,1,ikloc)) 
  wanc=zzero
  do n=1,nwantot
    do m=1,nwantot
      wanc(n,:)=wanc(n,:)+umtrx(m,n,ikloc)*wanncnrloc(m,:,ikloc)
    enddo
  enddo
  wanncnrloc(:,:,ikloc)=wanc(:,:)
  do n=1,nwantot
    do j=1,nstsv
      wann_unkmt(:,:,:,:,n,ikloc)=wann_unkmt(:,:,:,:,n,ikloc) + &
        &wfsvmtnrloc(:,:,:,:,j,ikloc)*wanncnrloc(n,j,ikloc)
      wann_unkit(:,:,n,ikloc)=wann_unkit(:,:,n,ikloc) + &
        &wfsvitnrloc(:,:,j,ikloc)*wanncnrloc(n,j,ikloc)
    enddo
  enddo
enddo
deallocate(wanc)
! calculate Wannier function occupancies 
wann_occ=0.d0
wann_ene=0.d0
do n=1,nwantot
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do j=1,nstsv
      w2=dreal(dconjg(wanncnrloc(n,j,ikloc))*wanncnrloc(n,j,ikloc))
      wann_occ(n)=wann_occ(n)+w2*occsvnr(j,ik)*wkptnr(ik)
      wann_ene(n)=wann_ene(n)+w2*evalsvnr(j,ik)*wkptnr(ik)
    enddo
  enddo
enddo
call mpi_grid_reduce(wann_occ(1),nwantot,dims=(/dim_k/),all=.true.)
call mpi_grid_reduce(wann_ene(1),nwantot,dims=(/dim_k/),all=.true.)
return
end subroutine


subroutine genwfnr(fout,lpmat)
use modmain
use mod_seceqn
implicit none
integer, intent(in) :: fout
logical, intent(in) :: lpmat
integer ik,ikloc,n,j,ik1,isym,i,ierr
complex(8), allocatable :: apwalm(:,:,:,:)
real(8) w2,t1,sz
logical, external :: wann_diel
complex(8), allocatable :: evecfvnrloc(:,:,:,:)
complex(8), allocatable :: evecsvnrloc(:,:,:)
complex(8), allocatable :: evecfdnrloc(:,:,:)
complex(8), allocatable :: evec(:,:)
!
!call gen_k_sym
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
  if (tsveqn) then
    sz=dble(nmatmax*nstfv*nspnfv)+dble(nstsv*nstsv)
  else
    sz=dble(nspinor*nmatmax*nstsv)
  endif
! wave-functions
  sz=sz+dble(lmmaxapw*nufrmax*natmtot*nstsv*nspinor)+dble(ngkmax*nstsv*nspinor)
! wannier functions
  if (wannier) then
    sz=sz+dble(nwantot*nstsv)+dble(lmmaxapw*nufrmax*natmtot*nspinor*nwantot)+&
      &dble(ngkmax*nspinor*nwantot)
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
allocate(wfsvmtnrloc(lmmaxapw,nufrmax,natmtot,nspinor,nstsv,nkptnrloc))
if (allocated(wfsvitnrloc)) deallocate(wfsvitnrloc)
allocate(wfsvitnrloc(ngkmax,nspinor,nstsv,nkptnrloc))
if (tsveqn) then
  allocate(evecfvnrloc(nmatmax,nstfv,nspnfv,nkptnrloc))
  allocate(evecsvnrloc(nstsv,nstsv,nkptnrloc))
else
  allocate(evecfdnrloc(nspinor*nmatmax,nstsv,nkptnrloc))
endif
if (lpmat) then
  if (allocated(pmatnrloc)) deallocate(pmatnrloc)
  allocate(pmatnrloc(3,nstsv,nstsv,nkptnrloc))
endif
if (wannier) then
  if (allocated(wanncnrloc)) deallocate(wanncnrloc)
  allocate(wanncnrloc(nwantot,nstsv,nkptnrloc))
  if (allocated(wann_unkmt)) deallocate(wann_unkmt)
  allocate(wann_unkmt(lmmaxapw,nufrmax,natmtot,nspinor,nwantot,nkptnrloc))
  if (allocated(wann_unkit)) deallocate(wann_unkit)
  allocate(wann_unkit(ngkmax,nspinor,nwantot,nkptnrloc))
endif
call timer_start(1,reset=.true.)
! read eigen-vectors
if (mpi_grid_side(dims=(/dim_k/))) then
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    if (tsveqn) then
      call getevecfv(vklnr(1,ik),vgklnr(1,1,ikloc),evecfvnrloc(1,1,1,ikloc))
      call getevecsv(vklnr(1,ik),evecsvnrloc(1,1,ikloc))
    else
      call getevecfd(vklnr(1,ik),vgklnr(1,1,ikloc),evecfdnrloc(1,1,ikloc)) 
    endif
  enddo !ikloc
endif !mpi_grid_side(dims=(/dim_k/)
call mpi_grid_barrier
! broadcast arrays
if (tsveqn) then
  call mpi_grid_bcast(evecfvnrloc(1,1,1,1),nmatmax*nstfv*nspnfv*nkptnrloc,&
    &dims=ortdims((/dim_k/)))
  call mpi_grid_bcast(evecsvnrloc(1,1,1),nstsv*nstsv*nkptnrloc,&
    &dims=ortdims((/dim_k/)))
else
  call mpi_grid_bcast(evecfdnrloc(1,1,1),nmatmax*nspinor*nstsv*nkptnrloc,&
    &dims=ortdims((/dim_k/)))
endif
call timer_stop(1)
if (wproc.and.fout.gt.0) then
  write(fout,'("Done in ",F8.2," seconds")')timer_get_value(1)
  if (fout.ne.6) call flushifc(fout)
endif
! generate wave functions from eigen vectors
wfsvmtnrloc=zzero
wfsvitnrloc=zzero
call timer_start(1,reset=.true.)
if (wproc.and.fout.gt.0) then
  write(fout,*)
  write(fout,'("Generating wave-functions")')
  if (fout.ne.6) call flushifc(fout)
endif
allocate(apwalm(ngkmax,lmmaxapw,apwordmax,natmtot))
allocate(evec(nspinor*nmatmax,nstsv))
if (ldisentangle.and..not.tsveqn) then
  write(*,*)
  write(*,'("Error(genwfnr): band disentanglment for the full diagonalization is not implemented")')
  call pstop
endif
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
! generate APW matching coefficients  
  call genapwalm(ngknr(ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc),&
      &sfacgknr(1,1,ikloc),apwalm)
  if (tsveqn) then
    call evecsvfd(evecfvnrloc(1,1,1,ikloc),evecsvnrloc(1,1,ikloc),evec)
  else
    evec(:,:)=evecfdnrloc(:,:,ikloc)
  endif
! generate wave-functions
  call genwfsvc(lmaxapw,lmmaxapw,ngknr(ikloc),nstsv,apwalm,&
    &evec,wfsvmtnrloc(1,1,1,1,1,ikloc),wfsvitnrloc(1,1,1,ikloc))
  if (wannier) then
    call wan_gencsv(lmmaxapw,vkcnr(1,ik),evalsvnr(1,ik),&
      &wfsvmtnrloc(1,1,1,1,1,ikloc),wanncnrloc(1,1,ikloc),ierr) 
    if (ierr.ne.0) then
      write(*,'("Warning(genwfnr): Wannier functions are wrong at k-point (ik, vkl) : ",I4,3G18.10)')ik,vklnr(:,ik)
    endif
    if (ldisentangle) then
! disentangle bands
      call disentangle(evalsvnr(1,ik),wanncnrloc(1,1,ikloc),&
        &evecsvnrloc(1,1,ikloc))
! generate wave-functions again
      call evecsvfd(evecfvnrloc(1,1,1,ikloc),evecsvnrloc(1,1,ikloc),evec)
      call genwfsvc(lmaxapw,lmmaxapw,ngknr(ikloc),nstsv,apwalm,&
        &evec,wfsvmtnrloc(1,1,1,1,1,ikloc),wfsvitnrloc(1,1,1,ikloc))
    endif
  endif
  if (lpmat) then
    call genpmatsv(ngknr(ikloc),igkignr(1,ikloc),vgkcnr(1,1,ikloc),&
      &wfsvmtnrloc(1,1,1,1,1,ikloc),wfsvitnrloc(1,1,1,ikloc),pmatnrloc(1,1,1,ikloc))
  endif
enddo !ikloc
deallocate(apwalm,evec)
call timer_stop(1)
if (wproc.and.fout.gt.0) then
  write(fout,'("Done in ",F8.2," seconds")')timer_get_value(1)
  call timestamp(fout)
  if (fout.ne.6) call flushifc(fout)
endif
! generate Bloch sums of Wannier functions
if (wannier) then
  call timer_start(1,reset=.true.)
  if (allocated(wann_unkmt)) deallocate(wann_unkmt)
  allocate(wann_unkmt(lmmaxapw,nufrmax,natmtot,nspinor,nwantot,nkptnrloc))
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
          &wfsvmtnrloc(:,:,:,:,j,ikloc)*wanncnrloc(n,j,ikloc)
        wann_unkit(:,:,n,ikloc)=wann_unkit(:,:,n,ikloc) + &
          &wfsvitnrloc(:,:,j,ikloc)*wanncnrloc(n,j,ikloc)
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
if (wproc.and.fout.gt.0) then
  write(fout,*)
  write(fout,'("band gap : ",F12.6)')bandgap
  write(fout,'("Ef       : ",F12.6)')efermi
  call timestamp(fout)
  if (fout.ne.6) call flushifc(fout)
endif
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
        wann_occ(n)=wann_occ(n)+w2*occsvnr(j,ik)*wkptnr(ik)
        wann_ene(n)=wann_ene(n)+w2*evalsvnr(j,ik)*wkptnr(ik)
      enddo
    enddo
  enddo
  call mpi_grid_reduce(wann_occ(1),nwantot,dims=(/dim_k/),all=.true.)
  call mpi_grid_reduce(wann_ene(1),nwantot,dims=(/dim_k/),all=.true.)
  if (wproc.and.fout.gt.0) then
    write(fout,'("  Wannier function occupancy and energy : ")')
    do n=1,nwantot
      write(fout,'("    n : ",I4,"  occ, ene : ",F8.6,2X,G18.10)')n,&
        &wann_occ(n),wann_ene(n)
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
if (spinpol.and.tsveqn) then
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
if (tsveqn) then
  deallocate(evecfvnrloc,evecsvnrloc)
else
  deallocate(evecfdnrloc)
endif
end subroutine

end module
