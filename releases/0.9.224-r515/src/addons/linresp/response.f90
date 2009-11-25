#ifdef _HDF5_
subroutine response
use modmain
#ifdef _MPI_
use mpi
#endif
use hdf5
implicit none

integer, allocatable :: ngknr(:)
integer, allocatable :: igkignr(:,:)
real(8), allocatable :: vgklnr(:,:,:)
real(8), allocatable :: vgkcnr(:,:,:)
real(8), allocatable :: gknr(:,:)
real(8), allocatable :: tpgknr(:,:,:)
complex(8), allocatable :: sfacgknr(:,:,:)

real(8), allocatable :: occsvnr(:,:)
real(8), allocatable :: evalsvnr(:,:)
complex(8), allocatable :: wfsvitloc(:,:,:,:)
complex(8), allocatable :: wfsvmtloc(:,:,:,:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfsvmt_t(:,:,:,:,:)
complex(8), allocatable :: wfsvit_t(:,:,:)
complex(8), allocatable :: wfc_t(:,:)
complex(8), allocatable :: pmat(:,:,:,:)

integer i,j,n,ik,ikloc,istsv,ik1,isym
integer sz,iint,iw
character*100 fname,qnm
character*3 c3
real(8) w2
integer, external :: iknrglob2
logical, external :: root_cart
logical, external :: in_cart

lmaxvr=4

! initialise universal variables
call init0
call init1
call pinit_cart

! MPI grid for tasks:
!  400 (matrix elements) : (1) k-points x (2) G-vectors x (3) q-points 
!  401 (chi0) : (1) k-points x (2) interband transition x (3) q-points 
!  402 (chi) : (1) energy mesh x (2) number of fxc kernels x (3) q-points

if (lrtype.eq.1.and..not.spinpol) then
  write(*,*)
  write(*,'("Error(response): can''t do magnetic response for unpolarized ground state")')
  write(*,*)
  call pstop
endif

if (.not.wannier) then
  lwannresp=.false.
  lwannopt=.false.
endif

if (task.eq.400) then
! read the density and potentials from file
  call readstate
! find the new linearisation energies
  call linengy
! generate the APW radial functions
  call genapwfr
! generate the local-orbital radial functions
  call genlofr
  call geturf
  call genurfprod
endif

if (task.eq.400.or.task.eq.401) then
  if (iproc.eq.0) call readfermi
  call dsync(efermi,1,.false.,.true.)
endif

if (in_cart()) then
  call qname(ivq0m_list(:,mpi_x(3)+1),qnm)
  if (root_cart((/1,1,0/))) call system("mkdir -p "//trim(qnm))
  call barrier(comm_cart)
  qnm="./"//trim(qnm)//"/"//trim(qnm)
  wproc=.false.
  if (task.eq.400) then
    if (root_cart((/1,1,0/))) then
      wproc=.true.
      fname=trim(qnm)//"_ME.OUT"
      open(150,file=trim(fname),form='formatted',status='replace')
    endif
  endif
  if (task.eq.401) then
    if (root_cart((/1,1,0/))) then
      wproc=.true.
      fname=trim(qnm)//"_CHI0.OUT"
      open(150,file=trim(fname),form='formatted',status='replace')
    endif
  endif
  if (task.eq.402) then
    write(c3,'(I3.3)')mpi_x(2)
    if (root_cart((/1,0,0/))) then
      wproc=.true.
      fname=trim(qnm)//"_CHI_"//c3//".OUT"
      open(150,file=trim(fname),form='formatted',status='replace')
    endif
  endif
  
  if (wproc) then
    write(150,'("Running on ",I8," proc.")')nproc
#ifdef _PIO_
    if (nproc.gt.1) then
      write(150,'("Using parallel I/O")')
    endif
#endif
    write(150,'("MPI grid size : ",3I6)')mpi_dims
    write(150,'("Wannier functions switch : ",L1)')wannier
    write(150,'("Response in local basis  : ",L1)')lwannresp
    call flushifc(150)
  endif
  
! distribute k-points over 1-st dimension of the grid
  if (allocated(nkptnrloc)) deallocate(nkptnrloc)
  allocate(nkptnrloc(0:mpi_dims(1)-1))
  if (allocated(ikptnrloc)) deallocate(ikptnrloc)
  allocate(ikptnrloc(0:mpi_dims(1)-1,2))
  call splitk(nkptnr,mpi_dims(1),nkptnrloc,ikptnrloc)
  nkptnr_loc=nkptnrloc(mpi_x(1))
  
   if (task.eq.400.or.task.eq.401.or.task.eq.404) then
! get occupancies and energies of states
    allocate(occsvnr(nstsv,nkptnr))
    allocate(evalsvnr(nstsv,nkptnr))
    call timer_start(3)
    if (wproc) then
      write(150,*)
      write(150,'("Reading energies and occupancies of states")')
      call flushifc(150)
    endif
    if (root_cart((/1,1,1/))) then
! read from IBZ
      do ik=1,nkpt
        call getoccsv(vkl(1,ik),occsv(1,ik))
        call getevalsv(vkl(1,ik),evalsv(1,ik))
      enddo
      do ik=1,nkptnr
        call findkpt(vklnr(1,ik),isym,ik1) 
        occsvnr(:,ik)=occsv(:,ik1)
        evalsvnr(:,ik)=evalsv(:,ik1)
      enddo
    endif
    call d_bcast_cart(comm_cart,evalsvnr,nstsv*nkptnr)
    call d_bcast_cart(comm_cart,occsvnr,nstsv*nkptnr)
    call timer_stop(3)
    if (wproc) then
      write(150,'("Done in ",F8.2," seconds")')timer(3,2)
      call flushifc(150)
    endif
  endif
  
  if (task.eq.400) then
! generate G+k vectors for entire BZ (this is required to compute 
!   wave-functions at each k-point)
    allocate(vgklnr(3,ngkmax,nkptnr_loc))
    allocate(vgkcnr(3,ngkmax,nkptnr_loc))
    allocate(gknr(ngkmax,nkptnr_loc))
    allocate(tpgknr(2,ngkmax,nkptnr_loc))
    allocate(ngknr(nkptnr_loc))
    allocate(sfacgknr(ngkmax,natmtot,nkptnr_loc))
    allocate(igkignr(ngkmax,nkptnr_loc))
    do ikloc=1,nkptnr_loc
      ik=iknrglob2(ikloc,mpi_x(1))
      call gengpvec(vklnr(1,ik),vkcnr(1,ik),ngknr(ikloc),igkignr(1,ikloc), &
        vgklnr(1,1,ikloc),vgkcnr(1,1,ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc))
      call gensfacgp(ngknr(ikloc),vgkcnr(1,1,ikloc),ngkmax,sfacgknr(1,1,ikloc))
    enddo
    allocate(wfsvmtloc(lmmaxvr,nrfmax,natmtot,nspinor,nstsv,nkptnr_loc))
    allocate(wfsvitloc(ngkmax,nspinor,nstsv,nkptnr_loc))
    allocate(evecfvloc(nmatmax,nstfv,nspnfv,nkptnr_loc))
    allocate(evecsvloc(nstsv,nstsv,nkptnr_loc))
    allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
    if (lwannopt) then
      allocate(pmat(3,nstsv,nstsv,nkptnr_loc))
    endif
    if (wproc) then
      sz=lmmaxvr*nrfmax*natmtot*nstsv*nspinor
      sz=sz+ngkmax*nstsv*nspinor
      sz=sz+nmatmax*nstfv*nspnfv
      sz=sz+nstsv*nstsv
      sz=16*sz*nkptnr_loc/1024/1024
      write(150,*)
      write(150,'("Size of wave-function arrays (MB) : ",I6)')sz
      write(150,*)
      write(150,'("Reading eigen-vectors")')
      call flushifc(150)
    endif
    call timer_start(1)
! read and transform eigen-vectors
    wfsvmtloc=zzero
    wfsvitloc=zzero
    if (root_cart((/0,1,1/))) then
#ifndef _PIO_
      do i=0,mpi_dims(1)-1
        if (i.eq.mpi_x(1)) then
#endif
          do ikloc=1,nkptnr_loc
            ik=iknrglob2(ikloc,mpi_x(1))
            call getevecfv(vklnr(1,ik),vgklnr(1,1,ikloc),evecfvloc(1,1,1,ikloc))
            call getevecsv(vklnr(1,ik),evecsvloc(1,1,ikloc))
! get apw coeffs 
            call match(ngknr(ikloc),gknr(1,ikloc),tpgknr(1,1,ikloc),        &
              sfacgknr(1,1,ikloc),apwalm)
! generate wave functions in muffin-tins
            call genwfsvmt(lmaxvr,lmmaxvr,ngknr(ikloc),evecfvloc(1,1,1,ikloc), &
              evecsvloc(1,1,ikloc),apwalm,wfsvmtloc(1,1,1,1,1,ikloc))
! generate wave functions in interstitial
            call genwfsvit(ngknr(ikloc),evecfvloc(1,1,1,ikloc), &
              evecsvloc(1,1,ikloc),wfsvitloc(1,1,1,ikloc))
            if (lwannopt) then
              call genpmat(ngknr(ikloc),igkignr(1,ikloc),vgkcnr(1,1,ikloc),&
                apwalm,evecfvloc(1,1,1,ikloc),evecsvloc(1,1,ikloc),pmat(1,1,1,ikloc))
            endif
          enddo !ikloc
#ifndef _PIO_
        endif
        call barrier(comm_cart_100)
      enddo
#endif
    endif !root_cart
    call barrier(comm_cart)
    call d_bcast_cart(comm_cart_011,wfsvmtloc, &
      lmmaxvr*nrfmax*natmtot*nspinor*nstsv*nkptnr_loc*2)
    call d_bcast_cart(comm_cart_011,wfsvitloc, &
      ngkmax*nspinor*nstsv*nkptnr_loc*2)
    call d_bcast_cart(comm_cart_011,evecfvloc, &
      nmatmax*nstfv*nspnfv*nkptnr_loc*2)
    call d_bcast_cart(comm_cart_011,evecsvloc, &
      nstsv*nstsv*nkptnr_loc*2)

! remove l content form particular bands      
!    if (.false.) then
!      allocate(wfsvmt_t(lmmaxvr,nrfmax,natmtot,nspinor,nstsv))
!      allocate(wfsvit_t(ngkmax,nspinor,nstsv))
!      do ikloc=1,nkptnr_loc
!        ik=iknrglob2(ikloc,mpi_x(1))
!        wfsvmt_t=wfsvmtloc(:,:,:,:,:,ikloc)
!        wfsvit_t=wfsvitloc(:,:,:,ikloc)
!        do j=1,nstsv
!          if (evalsvnr(j,ik).ge.-0.1d0.and.evalsvnr(j,ik).lt.0.1d0) then
!            wfsvmt_t(2:4,:,9:20,:,j)=0.d0
!          endif
!        enddo
!        wfsvmtloc(:,:,:,:,:,ikloc)=wfsvmt_t
!        wfsvitloc(:,:,:,ikloc)=wfsvit_t
!      enddo 
!      deallocate(wfsvmt_t)
!      deallocate(wfsvit_t)
!    endif
     if (.false.) then
       do ikloc=1,nkptnr_loc
         ik=iknrglob2(ikloc,mpi_x(1))
         do j=1,nstsv
           if (abs(evalsvnr(j,ik)-efermi).lt.0.01) then
             wfsvmtloc(:,:,:,:,j,ikloc)=zzero
             wfsvitloc(:,:,j,ikloc)=zzero
           endif
         enddo
       enddo
     endif  

    if (wannier) then
      if (allocated(wann_c)) deallocate(wann_c)
      allocate(wann_c(nwann,nstsv,nkptnr_loc))
      allocate(wfsvmt_t(lmmaxvr,nrfmax,natmtot,nspinor,nstsv))
      allocate(wfsvit_t(ngkmax,nspinor,nstsv))
      allocate(wfc_t(nwann,nstsv))
      if (wproc) then
        !write(150,*)
        write(150,'("  Generating Wannier functions")')
        if (lwfexpand) then
          write(150,*)
          if (laddwf) then
            write(150,'("    Adding Wannier functions content")')
          else
            write(150,'("    Removing Wannier functions content")')
            write(150,'("      dumping coefficient : ",F12.6)')alpha1
          endif
          do iint=1,nintwann      
            write(150,*)
            write(150,'("      energy interval : ",2F12.6)') &
              ewannint(1,iint),ewannint(2,iint)
            write(150,'("      Wannier functions : ",100I4)') &
              (iwannint(j,iint),j=1,nwannint(iint))
          enddo !iint          
        endif !lwfexpand
        call flushifc(150)
      endif !wproc
      do ikloc=1,nkptnr_loc
        ik=iknrglob2(ikloc,mpi_x(1))
        call genwann_c(ik,evalsvnr(1,ik),wfsvmtloc(1,1,1,1,1,ikloc),wfc_t)
        wann_c(:,:,ikloc)=wfc_t(:,:)
        if (lwfexpand) then
          if (laddwf) then
            wfsvmt_t=dcmplx(0.d0,0.d0)
            wfsvit_t=dcmplx(0.d0,0.d0)
          else
            wfsvmt_t=wfsvmtloc(:,:,:,:,:,ikloc)
            wfsvit_t=wfsvitloc(:,:,:,ikloc)
          endif
          do iint=1,nintwann      
            do j=1,nstsv
              if (evalsvnr(j,ik).ge.ewannint(1,iint).and.&
                  evalsvnr(j,ik).lt.ewannint(2,iint)) then
                do istsv=1,nstsv
                  do n=1,nwannint(iint)
                    iw=iwannint(n,iint)
                    if (laddwf) then
                      wfsvmt_t(:,:,:,:,j)=wfsvmt_t(:,:,:,:,j)+&
                        wfsvmtloc(:,:,:,:,istsv,ikloc)* &
                        wfc_t(iw,istsv)*dconjg(wfc_t(iw,j))
                      wfsvit_t(:,:,j)=wfsvit_t(:,:,j)+&
                        wfsvitloc(:,:,istsv,ikloc)*&
                        wfc_t(iw,istsv)*dconjg(wfc_t(iw,j))
                    else
                      wfsvmt_t(:,:,:,:,j)=wfsvmt_t(:,:,:,:,j)-&
                        wfsvmtloc(:,:,:,:,istsv,ikloc)* &
                        wfc_t(iw,istsv)*dconjg(wfc_t(iw,j))*alpha1
                      wfsvit_t(:,:,j)=wfsvit_t(:,:,j)-&
                        wfsvitloc(:,:,istsv,ikloc)*&
                        wfc_t(iw,istsv)*dconjg(wfc_t(iw,j))*alpha1
                    endif                       
                  enddo
                enddo
              endif
              enddo !j
          enddo !iint
          wfsvmtloc(:,:,:,:,:,ikloc)=wfsvmt_t
          wfsvitloc(:,:,:,ikloc)=wfsvit_t
        endif !lwfexpand
      enddo !ikloc
      deallocate(wfsvmt_t)
      deallocate(wfsvit_t)
      deallocate(wfc_t)
      wann_occ=0.d0
      do n=1,nwann
        do ikloc=1,nkptnr_loc
          ik=iknrglob2(ikloc,mpi_x(1))
          do j=1,nstsv
            w2=dreal(dconjg(wann_c(n,j,ikloc))*wann_c(n,j,ikloc))
            wann_occ(n)=wann_occ(n)+w2*occsvnr(j,ik)/nkptnr
          enddo
        enddo
      enddo
      call d_reduce_cart(comm_cart_100,.true.,wann_occ,nwann)
      wann_occ=wann_occ
      if (wproc) then
        write(150,'("    Wannier function occupation numbers : ")')
        do n=1,nwann
          write(150,'("      n : ",I4,"  occ : ",F8.6)')n,wann_occ(n)
        enddo
      endif
      lwanndiel=.true.
      do n=1,nwann
        if ((abs(wann_occ(n))*abs(wann_occ(n)-occmax)).gt.1d-8) &
          lwanndiel=.false.
      enddo
      if (wproc) then
        write(150,'("    Dielectric Wannier functions : ",L1)')lwanndiel
      endif
    endif !wannier
    
    call timer_stop(1)
    if (wproc) then
      write(150,'("Done in ",F8.2," seconds")')timer(1,2)
      call flushifc(150)
    endif
    deallocate(apwalm)
    deallocate(vgklnr)
    deallocate(vgkcnr)
    deallocate(gknr)
    deallocate(tpgknr)
    deallocate(sfacgknr)
  endif
  
  if (task.eq.400) then
! calculate matrix elements
    call timer_start(10)
    call response_me(ivq0m_list(1,mpi_x(3)+1),wfsvmtloc,wfsvitloc,ngknr, &
      igkignr,occsvnr,evalsvnr,pmat)
    call timer_stop(10)
    if (wproc) then
      write(150,'("Total time : ",F8.2," seconds")')timer(10,2)
      call flushifc(150)
    endif
  endif
  
  if (task.eq.401) then
! calculate chi0
    call timer_start(11)
    call response_chi0(ivq0m_list(1,mpi_x(3)+1),evalsvnr)
    call timer_stop(11)
    if (wproc) then
      write(150,'("Total time : ",F8.2," seconds")')timer(11,2)
      call flushifc(150)
    endif
  endif
  
  if (task.eq.402) then
! calculate chi
    call timer_start(12)
    call response_chi(ivq0m_list(1,mpi_x(3)+1))
    call timer_stop(12)
    if (wproc) then
      write(150,'("Total time : ",F8.2," seconds")')timer(12,2)
      call flushifc(150)
    endif
  endif
  
  if (task.eq.404) call response_jdos(occsvnr,evalsvnr)
  
  if (wproc) close(150)
  
  if (task.eq.400.and.lwannopt) deallocate(pmat)
  
  if (task.eq.400.or.task.eq.403) then
    deallocate(wfsvmtloc)
    deallocate(wfsvitloc)
    deallocate(ngknr)
    deallocate(igkignr)
    deallocate(occsvnr)
  endif

endif !in_cart()

return
end
#endif