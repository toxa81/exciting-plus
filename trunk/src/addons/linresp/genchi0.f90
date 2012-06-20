subroutine genchi0(iq)
use modmain
use mod_addons_q
use mod_nrkp
use mod_wannier
use mod_expigqr
use mod_linresp
implicit none
! arguments
integer, intent(in) :: iq
! local variables
complex(8), allocatable :: chi0(:,:)
complex(8), allocatable :: chi0wan_k(:,:,:)
complex(8), allocatable :: chi0wan(:,:)
complex(8), allocatable :: mexp(:,:,:)
complex(8), allocatable :: megqwan_tmp(:,:)
integer i,iw,i1,i2,ikloc,n,j
integer ist1,ist2,nwloc,jwloc
integer it1(3),it2(3),it(3)
character*100 qnm,qdir,fout,fstat
integer ik,n1,n2
integer nwt
integer, allocatable :: iwt_tmp(:,:)
real(8) vtc(3)

real(8) t1,t2,t3

call getqdir(iq,vqm(:,iq),qdir)
call getqname(vqm(:,iq),qnm)
qnm=trim(qdir)//"/"//trim(qnm)
wproc=.false.
if (mpi_grid_root((/dim_k/))) then
  wproc=.true.
  fout=trim(qnm)//"_CHI0.OUT"
  open(150,file=trim(fout),form="FORMATTED",status="REPLACE")
  fstat=trim(qnm)//"_chi0_stat.txt"
endif

call papi_timer_start(pt_chi0)

if (wproc) then
  write(150,*)
  write(150,'("Calculation of chi0")')
  write(150,*)
  write(150,'("Energy mesh parameters:")')
  if (timgw) then
    write(150,'("  number of imaginary frequencies : ",I4)')lr_nw
    write(150,'("  frequency interval [eV] : ", 2F9.4)')lr_iw0,lr_iw1
  else
    write(150,'("  energy interval [eV] : ", 2F9.4)')lr_w0*ha2ev,lr_w1*ha2ev
    write(150,'("  energy step     [eV] : ", F9.4)')lr_dw*ha2ev
    write(150,'("  eta             [eV] : ", F9.4)')lr_eta*ha2ev
  endif
  write(150,*)  
  write(150,'("Included band interval (Ha)        : ",2F8.2)')&
    &chi0_include_bands(1),chi0_include_bands(2)
  write(150,'("Excluded band interval (Ha)        : ",2F8.2)')&
    &chi0_exclude_bands(1),chi0_exclude_bands(2) 
  call flushifc(150)
endif
  
if (wannier_chi0_chi) then
! Warning: this is quick fix. Code is taken from genmegq.f90.  
!   The reason for recomputing megqwan is that for optical limit of 
!   response in Wannier basis we need matrix elements with correct 
!   asymptotics for q->0. This was done for megqblh in genmegq.f90,
!   but megqwan was computed before fixing the q->0 limit (for cRPA 
!   calculations we need those "normal" megqwan elements). 
! There is still something wrong with optical limit of Wannier 
!   response for metallic systems
  if (all(vqm(:,iq).eq.0)) then
    megqwan=zzero
    call genmegqwan(iq)
    call mpi_grid_reduce(megqwan(1,1),megqwantran%nwt*ngq(iq),&
      &dims=(/dim_k/),all=.true.)
    megqwan=megqwan/nkptnr
  endif
endif
! filter matrix elements
if (wannier_chi0_chi) then
  allocate(megqwan_tmp(megqwantran%nwt,ngq(iq)))
  megqwan_tmp=zzero
  allocate(iwt_tmp(5,megqwantran%nwt))
  iwt_tmp=0
  nwt=0
  do i=1,megqwantran%nwt
    if (abs(megqwan(i,iig0q)).ge.megqwan_cutoff(1).and.&
        &abs(megqwan(i,iig0q)).le.megqwan_cutoff(2)) then
      nwt=nwt+1
      iwt_tmp(:,nwt)=megqwantran%iwt(:,i)
      megqwan_tmp(nwt,:)=megqwan(i,:)
    endif
  enddo
  call deletewantran(megqwantran)
  megqwantran%nwt=nwt
  allocate(megqwantran%iwt(5,nwt))
  megqwantran%iwt(:,1:nwt)=iwt_tmp(:,1:nwt)
  deallocate(megqwan)
  allocate(megqwan(nwt,ngq(iq)))
  megqwan(1:nwt,:)=megqwan_tmp(1:nwt,:)
  deallocate(iwt_tmp,megqwan_tmp)
endif
!if (wannier_chi0_chi.and..true.) then
!  do i=1,megqwantran%nwt 
!    if (.not.(megqwantran%iwt(1,i).le.5.and.&
!              megqwantran%iwt(2,i).le.5.and.&
!              all(megqwantran%iwt(3:5,i).eq.0))) megqwan(i,:)=zzero
!  enddo
!endif

! for response in Wannier bais
if (wannier_chi0_chi) then
  if (wproc) then
    write(150,*)
    write(150,'("Wannier chi0 AFM : ",L1)')wannier_chi0_afm
    write(150,*)
    write(150,'("megqwan value cutoff (min,max) : ",2F12.6)')&
      &megqwan_cutoff
    write(150,'("megqwan distance cutoff (min,max) : ",2F12.6)')&
      &megqwan_mindist,megqwan_maxdist
    write(150,'("Number of Wannier transitions after cutoff : ",I6)')&
      &megqwantran%nwt
    write(150,*)
    write(150,'("<n|e^{-iqx}|n''T> where q-vector is not reduced")')    
    write(150,'(65("-"))')    
    call printwanntrans(150,megqwan(1,iig0q))
  endif
  allocate(chi0wan(megqwantran%nwt,megqwantran%nwt))
  allocate(chi0wan_k(megqwantran%nwt,megqwantran%nwt,nkptnrloc))
  allocate(mexp(megqwantran%nwt,megqwantran%nwt,nkptnrloc))
  do i1=1,megqwantran%nwt
    do i2=1,megqwantran%nwt
      it1(:)=megqwantran%iwt(3:5,i1)
      it2(:)=megqwantran%iwt(3:5,i2)
      it(:)=it1(:)-it2(:)
      vtc(:)=avec(:,1)*it(1)+avec(:,2)*it(2)+avec(:,3)*it(3)
      do ikloc=1,nkptnrloc
        ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
! phase e^{i(k+q)(T1-T2)}
        mexp(i1,i2,ikloc)=exp(zi*dot_product(vkcnr(:,ik)+vqc(:,iq),vtc(:)))
      enddo
    enddo
  enddo
! arrangement for zgemm  
  allocate(wann_cc(nmegqblhwanmax,megqwantran%nwt,nkptnrloc))
  allocate(wann_cc2(nmegqblhwanmax,megqwantran%nwt))
  wann_cc=zzero
  do ikloc=1,nkptnrloc
    do i1=1,nmegqblhwan(ikloc)
      i=imegqblhwan(i1,ikloc)
      ist1=bmegqblh(1,i,ikloc)
      ist2=bmegqblh(2,i,ikloc)
      do n=1,megqwantran%nwt
        n1=megqwantran%iwt(1,n)
        n2=megqwantran%iwt(2,n)
        wann_cc(i1,n,ikloc)=wanncnrloc(n1,ist1,ikloc)*&
          &dconjg(wann_c_jk(n2,ist2,ikloc))
      enddo
    enddo !i1
  enddo !ikloc
endif !wannier_chi0_chi

allocate(chi0(ngq(iq),ngq(iq)))
allocate(megqblh2(nstsv*nstsv,ngq(iq)))

! distribute frequency points over 1-st dimension
nwloc=mpi_grid_map(lr_nw,dim_k)
if (allocated(chi0loc)) deallocate(chi0loc)
allocate(chi0loc(ngq(iq),ngq(iq),nwloc))
if (wannier_chi0_chi) then
  if (allocated(chi0wanloc)) deallocate(chi0wanloc)
  allocate(chi0wanloc(megqwantran%nwt,megqwantran%nwt,nwloc))
endif

call timer_start(1,reset=.true.)
call timer_reset(2)
call timer_reset(3)
! loop over energy points
do iw=1,lr_nw
  call timer_start(10,reset=.true.)
! sum over fraction of k-points
  call timer_start(2)
  chi0=zzero
  do ikloc=1,nkptnrloc
    if (nmegqblh(ikloc).gt.0) then
! for each k-point : sum over interband transitions
      call genchi0blh(ikloc,ngq(iq),lr_w(iw),chi0)
    endif
  enddo
! find the processor j which will get the full chi0 and chi0wan matrices
  jwloc=mpi_grid_map(lr_nw,dim_k,glob=iw,x=j)
! sum over k-points and band transitions
  call mpi_grid_reduce(chi0(1,1),ngq(iq)*ngq(iq),dims=(/dim_k/),&
    &root=(/j/))
  chi0=chi0/nkptnr/omega
! processor j saves chi0 to local array  
! TODO: with new API this can be simplified
  if (mpi_grid_dim_pos(dim_k).eq.j) chi0loc(:,:,jwloc)=chi0(:,:)
  call timer_stop(2)
! for response in Wannier basis
  if (wannier_chi0_chi) then
    call timer_start(3)
    chi0wan_k=zzero
    chi0wan=zzero
    do ikloc=1,nkptnrloc
      if (nmegqblhwan(ikloc).gt.0) then
        call genchi0wan_k(ikloc,lr_w(iw),chi0wan_k(1,1,ikloc))
      endif
      chi0wan(:,:)=chi0wan(:,:)+mexp(:,:,ikloc)*chi0wan_k(:,:,ikloc)
    enddo !ikloc
! sum chi0wan over all k-points
    call mpi_grid_reduce(chi0wan(1,1),megqwantran%nwt*megqwantran%nwt,&
      &dims=(/dim_k/),root=(/j/))
    chi0wan(:,:)=chi0wan(:,:)/nkptnr/omega
    if (wannier_chi0_afm) chi0wan(:,:)=chi0wan(:,:)*2.d0
!!! WSTHORNTON
    !call svdchi0(chi0wan,iq,iw,lr_w(iw)*ha2ev)
! processor j saves chi0wan to local array  
    if (mpi_grid_dim_pos(dim_k).eq.j) chi0wanloc(:,:,jwloc)=chi0wan(:,:)
    call timer_stop(3)
  endif !wannier_chi0_chi
  call timer_stop(10)
  if (wproc) then
    open(160,file=trim(fstat),status="REPLACE",form="FORMATTED")
    write(160,'(I8)')iw
    write(160,'(F8.2)')timer_get_value(10)   
    close(160)
  endif
enddo !iw
call timer_stop(1)
t1=timer_get_value(1)
t2=timer_get_value(2)
t3=timer_get_value(3)
if (wproc) then
  write(150,*)
  write(150,'("Total time per frequency point   : ",F8.2)')t1/lr_nw
  write(150,'("  Bloch basis part (chi0)        : ",F8.2)')t2/lr_nw
  write(150,'("  Wannier basis part (chi0)      : ",F8.2)')t3/lr_nw
  call flushifc(150)
endif

call papi_timer_stop(pt_chi0)

call mpi_grid_barrier(dims=(/dim_k/))

deallocate(chi0)
deallocate(megqblh2)
if (wannier_chi0_chi) then
  deallocate(chi0wan)
  deallocate(chi0wan_k)
  deallocate(mexp)
  deallocate(wann_cc)
  deallocate(wann_cc2)
endif
if (wproc) then
  write(150,*)
  write(150,'("Done.")')
  call flushifc(150)
  close(150)
endif
return
end
