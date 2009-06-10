subroutine dosrlm
use modmain
implicit none
! local variables
integer lmax,lmmax,l,m,lm,nsk(3)
integer ik,ispn,is,ia,ias,ist,iw,i,ikloc,idx0,bs
real(8) t1
! allocatable arrays
real(4), allocatable :: bndchr(:,:,:,:,:)
!real(8), allocatable :: e(:,:,:)
real(8), allocatable :: f(:,:)
real(8), allocatable :: w(:)
real(8), allocatable :: tdos(:,:)
real(8), allocatable :: idos(:,:)
real(8), allocatable :: pdos(:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: sdmat(:,:,:,:)
integer, external :: ikglob

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

call geturf
call genurfprod

lmax=min(3,lmaxapw)
lmmax=(lmax+1)**2
allocate(bndchr(lmmax,natmtot,nspinor,nstsv,nkpt))
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
allocate(sdmat(nspinor,nspinor,nstsv,nkpt))

evalsv=0.d0
bndchr=0.0
sdmat=dcmplx(0.d0,0.d0)
do i=0,nproc-1
  if (iproc.eq.i) then
    do ikloc=1,nkptloc(iproc)
      ik=ikglob(ikloc)
      call getevecfv(vkl(1,ik),vgkl(1,1,1,ikloc),evecfv)
      call getevecsv(vkl(1,ik),evecsv)
      call getevalsv(vkl(1,ik),evalsv(1,ik))
! substract Fermi energy
      evalsv(:,ik)=evalsv(:,ik)-efermi
! correction for scissors operator
      do ist=1,nstsv
        if (evalsv(ist,ik).gt.0.d0) evalsv(ist,ik)=evalsv(ist,ik)+scissor
      enddo
! compute the band character
      call bandchar(.true.,lmax,ikloc,evecfv,evecsv,lmmax,bndchr(1,1,1,1,ik))
!      if (wannier) then
!        call getwfc(ik,wann_c(1,1,1,ik)) 
!      endif
      call gensdmat(evecsv,sdmat(:,:,:,ik))
    end do
  endif
  call barrier(comm_world)
enddo
do ik=1,nkpt
  call rsync(bndchr(1,1,1,1,ik),lmmax*natmtot*nspinor*nstsv,.true.,.true.)
  call zsync(sdmat(1,1,1,ik),nspinor*nspinor*nstsv,.true.,.true.)
enddo
call dsync(evalsv,nstsv*nkpt,.true.,.true.)
  
! generate energy grid
wdos(1)=minval(evalsv(:,:))-0.1
wdos(2)=maxval(evalsv(:,:))+0.1
t1=0.001d0
nwdos=1+(wdos(2)-wdos(1))/t1
allocate(w(nwdos))
do iw=1,nwdos
  w(iw)=t1*dble(iw-1)+wdos(1)
end do

! allocate local arrays
allocate(f(nstsv,nkpt))
allocate(pdos(nwdos,lmmax,nspinor,natmtot))
if (iproc.eq.0) then
  allocate(tdos(nwdos,nspinor))
  allocate(idos(nwdos,nspinor))
endif
! number of subdivisions used for interpolation
do i=1,3
  nsk(i)=max(ngrdos/ngridk(i),1)
end do
! compute total DOS
if (iproc.eq.0) then
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
call idxbos(natmtot,nproc,iproc+1,idx0,bs)
do ias=idx0+1,idx0+bs
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
  call dsync(pdos(1,1,1,ias),nwdos*lmmax*nspinor,.true.,.false.)
enddo

if (iproc.eq.0) then
! write total DOS
  open(50,file='TDOS_EV.OUT',action='WRITE',form='FORMATTED')
  do ispn=1,nspinor
    if (ispn.eq.1) then
      t1=1.d0
    else
      t1=-1.d0
    end if
    do iw=1,nwdos
      write(50,'(2G18.10)') w(iw)*ha2ev,t1*tdos(iw,ispn)/ha2ev
    end do
    write(50,'("     ")')
  enddo
  close(50)
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
  open(50,file='IDOS_EV.OUT',action='WRITE',form='FORMATTED')
  do ispn=1,nspinor
    if (ispn.eq.1) then
      t1=1.d0
    else
      t1=-1.d0
    end if
    do iw=1,nwdos
      write(50,'(2G18.10)') w(iw)*ha2ev,t1*idos(iw,ispn)/ha2ev
    end do
    write(50,'("     ")')  
  end do
  close(50)
! write partial DOS
  open(50,file='PDOS.OUT',form='UNFORMATTED',status='REPLACE')
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
    write(50)w(iw)*ha2ev
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
  close(50)
endif
deallocate(f,w,pdos)
if (iproc.eq.0) then
  deallocate(tdos,idos)
endif
!if (wannier) then
!  do n=1,wf_dim
!    do ik=1,nkpt
!      do ist=1,nstfv
!        f(ist,ik)=abs(wfc(n,ist,1,ik))**2
!      end do
!    end do
!    call brzint(nsmdos,ngridk,nsk,ikmap,nwdos,wdos,nstsv,nstsv, &
!      e(1,1,1),f,gp) 
!    write(fname,'("WDOS_",I4.4,".OUT")')n
!    open(50,file=trim(fname),action='WRITE',form='FORMATTED')
!    do iw = 1, nwdos
!      write(50,'(2G18.10)')w(iw)*ha2ev, gp(iw)/ha2ev
!    enddo
!    close(50)
!  enddo
!endif 
!endif
deallocate(bndchr,evecfv,evecsv)
return
end subroutine
!EOC
