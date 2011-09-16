subroutine dosrlm
use modmain
use mod_wannier
use mod_seceqn
implicit none
! local variables
integer lmax,lmmax,l,m,lm,nsk(3)
integer ik,ispn,is,ia,ias,ist,iw,i,ikloc,x0,n
real(8) t1
complex(8) zt1
! allocatable arrays
real(4), allocatable :: bndchr(:,:,:,:,:)
real(8), allocatable :: f(:,:)
real(8), allocatable :: w(:)
real(8), allocatable :: tdos(:,:)
real(8), allocatable :: idos(:,:)
real(8), allocatable :: pdos(:,:,:,:)
real(8), allocatable :: doswan(:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: evecfd(:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfsvmt(:,:,:,:,:) 
complex(8), allocatable :: sdmat(:,:,:,:)
integer iasloc,natmtotloc

! initialise universal variables
call init0
call init1

! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr
! get radial-muffint tin functions
call getufr
! get product of radial functions
call genufrp  

lmax=min(3,lmaxapw)
lmmax=(lmax+1)**2
allocate(bndchr(lmmax,natmtot,nspinor,nstsv,nkpt))
allocate(sdmat(nspinor,nspinor,nstsv,nkpt))
allocate(wfsvmt(lmmax,nufrmax,natmtot,nspinor,nstsv))
if (tsveqn) then
  allocate(evecfv(nmatmax,nstfv,nspnfv))
  allocate(evecsv(nstsv,nstsv))
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
else
  allocate(evecfd(nspinor*nmatmax,nstsv))
  allocate(apwalm(ngkmax,lmmaxapw,apwordmax,natmtot))
endif
evalsv=0.d0
bndchr=0.0
sdmat=zzero
if (mpi_grid_side(dims=(/dim_k/))) then
  do ikloc=1,nkptloc
    ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
    call getevalsv(vkl(1,ik),evalsv(1,ik))
    if (tsveqn) then
      call getevecfv(vkl(1,ik),vgkl(1,1,1,ikloc),evecfv)
      call getevecsv(vkl(1,ik),evecsv)
      call match(ngk(1,ik),gkc(1,1,ikloc),tpgkc(1,1,1,ikloc),&
        sfacgk(1,1,1,ikloc),apwalm)
      call genwfsvmt(lmax,lmmax,ngk(1,ik),evecfv(:,:,1),evecsv,apwalm,wfsvmt)
      if (wannier) call genwann(ikloc,evecfv,evecsv)
      call gensdmat(evecsv,sdmat(:,:,:,ik))
    else
      call getevecfd(vkl(1,ik),vgkl(1,1,1,ikloc),evecfd)
      call genapwalm(ngk(1,ik),gkc(1,1,ikloc),tpgkc(1,1,1,ikloc),&
        sfacgk(1,1,1,ikloc),apwalm)
      call genwfsvc(lmax,lmmax,ngk(1,ik),nstsv,apwalm,evecfd,wfsvmt)
      if (wannier) call wan_gencsv(lmmax,vkc(1,ik),evalsv(1,ik),&
        wfsvmt,wann_c(1,1,ikloc))
      call gensdmatfd(nmat(1,ik),evecfd,sdmat(:,:,:,ik))
    endif
! compute the band character
    call bandchar(.true.,lmax,lmmax,wfsvmt,bndchr(1,1,1,1,ik))
  end do
endif
do ik=1,nkpt
  call mpi_grid_reduce(bndchr(1,1,1,1,ik),lmmax*natmtot*nspinor*nstsv,&
    dims=(/dim_k/),side=.true.,all=.true.)
  call mpi_grid_reduce(sdmat(1,1,1,ik),nspinor*nspinor*nstsv,&
    dims=(/dim_k/),side=.true.,all=.true.)
enddo
call mpi_grid_reduce(evalsv(1,1),nstsv*nkpt,dims=(/dim_k/),side=.true.,&
  all=.true.)
  
! generate energy grid
wdos(1)=minval(evalsv(:,:))-0.1
wdos(2)=maxval(evalsv(:,:))+0.1
wdos(1)=int(wdos(1)*1.0d3)*1.0d-3
wdos(2)=int(wdos(2)*1.0d3)*1.0d-3
t1=1.0d-3
nwdos=1+(wdos(2)-wdos(1))/t1
allocate(w(nwdos))
do iw=1,nwdos
  w(iw)=t1*dble(iw-1)+wdos(1)
end do

! allocate local arrays
allocate(f(nstsv,nkpt))
allocate(pdos(nwdos,lmmax,nspinor,natmtot))
if (mpi_grid_root()) then
  allocate(tdos(nwdos,nspinor))
  allocate(idos(nwdos,nspinor))
endif
! number of subdivisions used for interpolation
do i=1,3
  nsk(i)=max(ngrdos/ngridk(i),1)
end do
! compute total DOS
if (mpi_grid_root()) then
  do ispn=1,nspinor
    do ik=1,nkpt
      do ist=1,nstsv
! use spin character for weight
        f(ist,ik)=dble(sdmat(ispn,ispn,ist,ik))
      end do
    end do !ik
    call brzint(nsmdos,ngridk,nsk,ikmap,nwdos,wdos,nstsv,nstsv,evalsv,f, &
      tdos(1,ispn))
  end do !ispn
  tdos=tdos*occmax
endif
! compute partial DOS in parallel
pdos=0.d0
natmtotloc=mpi_grid_map(natmtot,dim_k)
do iasloc=1,natmtotloc
  ias=mpi_grid_map(natmtot,dim_k,loc=iasloc)
  do ispn=1,nspinor
    do l=0,lmax
      do m=-l,l
        lm=idxlm(l,m)
        do ik=1,nkpt
          do ist=1,nstsv
            f(ist,ik)=bndchr(lm,ias,ispn,ist,ik)
          end do
        end do
        call brzint(nsmdos,ngridk,nsk,ikmap,nwdos,wdos,nstsv,nstsv, &
         evalsv,f,pdos(1,lm,ispn,ias))
      end do !m
    end do !l
  end do !ispn
enddo !ias
pdos=pdos*occmax
do ias=1,natmtot
  call mpi_grid_reduce(pdos(1,1,1,ias),nwdos*lmmax*nspinor,dims=(/dim_k/),&
    side=.true.)
enddo
allocate(doswan(nwdos,nwantot))
doswan=0.d0
if (wannier) then
  do n=1,nwantot
    f=0.d0
    do ik=1,nkpt
      ikloc=mpi_grid_map(nkpt,dim_k,x=x0,glob=ik)
      if (mpi_grid_dim_pos(dim_k).eq.x0) then
        do ist=1,nstsv
          f(ist,ik)=abs(wann_c(n,ist,ikloc))**2
        end do
      endif
    end do
    call mpi_grid_reduce(f(1,1),nstsv*nkpt,dims=(/dim_k/))
    call brzint(nsmdos,ngridk,nsk,ikmap,nwdos,wdos,nstsv,nstsv, &
      evalsv,f,doswan(1,n)) 
  enddo
endif 

if (mpi_grid_root()) then
! write total DOS
  open(50,file="TDOS.OUT",action="WRITE",form="FORMATTED")
  open(51,file="tdos.dat",action="WRITE",form="FORMATTED")
  do ispn=1,nspinor
    if (ispn.eq.1) then
      t1=1.d0
    else
      t1=-1.d0
    end if
    do iw=1,nwdos
      write(50,'(2G18.10)') w(iw),t1*tdos(iw,ispn)
      write(51,'(2G18.10)') (w(iw)-efermi)*ha2ev,t1*tdos(iw,ispn)/ha2ev
    end do
    write(50,'("     ")')
    write(51,'("     ")')
  enddo
  close(50)
  close(51)
! compute interstitial DOS
  idos=tdos
  do ispn=1,nspinor
    do ias=1,natmtot
      do lm=1,lmmax
        idos(:,ispn)=idos(:,ispn)-pdos(:,lm,ispn,ias)
      enddo
    enddo
  enddo
! write interstitial DOS
  open(50,file="IDOS.OUT",action="WRITE",form="FORMATTED")
  open(51,file="idos.dat",action="WRITE",form="FORMATTED")
  do ispn=1,nspinor
    if (ispn.eq.1) then
      t1=1.d0
    else
      t1=-1.d0
    end if
    do iw=1,nwdos
      write(50,'(2G18.10)') w(iw),t1*idos(iw,ispn)
      write(51,'(2G18.10)') (w(iw)-efermi)*ha2ev,t1*idos(iw,ispn)/ha2ev
    end do
    write(50,'("     ")')  
    write(51,'("     ")')  
  end do
  close(50)
  close(51)
! write partial DOS
  open(50,file="PDOS.OUT",form="UNFORMATTED",status="REPLACE")
  write(50)nspecies,natmtot,nspinor,lmax,lmmax,nwdos
  do is=1,nspecies
    write(50)spsymb(is)
    write(50)natoms(is)
  enddo
  do is=1,nspecies
    do ia=1,natoms(is)
      write(50)idxas(ia,is),is
    enddo
  enddo
  do iw=1,nwdos
    write(50)(w(iw)-efermi)*ha2ev
  enddo
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ispn=1,nspinor
        do l=0,lmax
          do m=-l,l
            lm=idxlm(l,m)
            do iw=1,nwdos
              write(50)pdos(iw,lm,ispn,ias)/ha2ev
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
  write(50)wannier
  if (wannier) then
    write(50)nwantot
    do n=1,nwantot
      do iw=1,nwdos
        write(50)occmax*doswan(iw,n)/ha2ev
      enddo
    enddo
  endif
  close(50)
endif
if (mpi_grid_root()) then
  deallocate(tdos,idos)
endif
deallocate(f,w,pdos,doswan)
if (tsveqn) then
  deallocate(evecfv,evecsv)
else
  deallocate(evecfd)
endif
deallocate(bndchr,apwalm,wfsvmt)
return
end subroutine
!EOC
