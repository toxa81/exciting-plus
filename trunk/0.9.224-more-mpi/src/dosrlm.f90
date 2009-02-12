subroutine dosrlm
use modmain
implicit none
! local variables
integer lmax,lmmax,l,m,lm,nsk(3)
integer ik,ispn,is,ia,ias,ist,iw,i
real(8) t1
! allocatable arrays
real(4), allocatable :: bndchr(:,:,:,:,:)
real(8), allocatable :: e(:,:,:)
real(8), allocatable :: f(:,:)
real(8), allocatable :: w(:)
real(8), allocatable :: g(:,:)
real(8), allocatable :: gp(:)
real(8), allocatable :: pdos(:,:)
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: sdmat(:,:,:,:)


! initialise universal variables
call init0
call init1

lmax=min(3,lmaxapw)
lmmax=(lmax+1)**2

allocate(bndchr(lmmax,natmtot,nspinor,nstsv,nkpt))
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
allocate(sdmat(nspinor,nspinor,nstsv,nkpt))

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

do ik=1,nkpt
  call getevecfv(vkl(1,ik),vgkl(1,1,1,ik),evecfv)
  call getevecsv(vkl(1,ik),evecsv)
  call getevalsv(vkl(1,ik),evalsv(1,ik))
  call bandchar(.true.,lmax,ik,evecfv,evecsv,lmmax,bndchr(1,1,1,1,ik))
  if (wannier) then
    call getwfc(ik,wann_c(1,1,1,ik)) 
  endif
  call gensdmat(evecsv,sdmat(:,:,:,ik))
end do

! generate energy grid
wdos(1)=minval(evalsv(:,:)-efermi)-0.1
wdos(2)=maxval(evalsv(:,:)-efermi)+0.1
t1=0.001d0
!t1=(wdos(2)-wdos(1))/dble(nwdos)
nwdos=1+(wdos(2)-wdos(1))/t1

! allocate local arrays
allocate(e(nstsv,nkpt,nspinor))
allocate(f(nstsv,nkpt))
allocate(w(nwdos))
allocate(g(nwdos,nspinor))
allocate(gp(nwdos))
allocate(pdos(nwdos,lmmax))

do iw=1,nwdos
  w(iw)=t1*dble(iw-1)+wdos(1)
end do
! number of subdivisions used for interpolation
do i=1,3
  nsk(i)=max(ngrdos/ngridk(i),1)
end do

! write total DOS
open(50,file='TDOS_EV.OUT',action='WRITE',form='FORMATTED')
do ispn=1,nspinor
  if (ispn.eq.1) then
    t1=1.d0
  else
    t1=-1.d0
  end if
  do ik=1,nkpt
    do ist=1,nstsv
! subtract the Fermi energy
      e(ist,ik,ispn)=evalsv(ist,ik)-efermi
! correction for scissors operator
      if (e(ist,ik,ispn).gt.0.d0) e(ist,ik,ispn)=e(ist,ik,ispn)+scissor
! use spin character for weight
      f(ist,ik)=dble(sdmat(ispn,ispn,ist,ik))
    end do
  end do !ik
  call brzint(nsmdos,ngridk,nsk,ikmap,nwdos,wdos,nstsv,nstsv,e(1,1,ispn),f, &
   g(1,ispn))
  g=occmax*g 
  do iw=1,nwdos
    write(50,'(2G18.10)') w(iw)*ha2ev,t1*g(iw,ispn)/ha2ev
  end do
  write(50,'("     ")')
end do !ispn
close(50)

! write partial DOS
open(160,file='PDOS.OUT',form='UNFORMATTED',status='REPLACE')
write(160)nspecies,natmtot,nspinor,lmax,lmmax,nwdos
do is=1,nspecies
  write(160)spsymb(is)
  write(160)natoms(is)
enddo
do is=1,nspecies
  do ia=1,natoms(is)
    write(160)idxas(ia,is),is
  enddo
enddo
do iw=1,nwdos
  write(160)w(iw)*ha2ev
enddo
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
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
           e(1,1,ispn),f,gp)
          gp=occmax*gp 
          do iw=1,nwdos
            write(160)gp(iw)/ha2ev
            g(iw,ispn)=g(iw,ispn)-gp(iw)
          end do            
        end do !m
      end do !l
    end do !ispn
  end do !ia
end do !is
close(160)

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

! write interstitial DOS
open(50,file='IDOS_EV.OUT',action='WRITE',form='FORMATTED')
do ispn=1,nspinor
  if (ispn.eq.1) then
    t1=1.d0
  else
    t1=-1.d0
  end if
  do iw=1,nwdos
    write(50,'(2G18.10)') w(iw)*ha2ev,t1*g(iw,ispn)/ha2ev
  end do
  write(50,'("     ")')  
end do
close(50)

deallocate(bndchr,e,f,w,g,gp,sdmat)
deallocate(evecfv,evecsv)
return
end subroutine
!EOC
