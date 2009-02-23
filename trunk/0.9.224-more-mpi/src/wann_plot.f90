subroutine wann_plot
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none

real(8) r(3),t(3),vr0l(3)
real(8) bound3d(3,3),orig3d(3),zero3d(3)
real(8) bound2d(3,2),orig2d(3)
complex(4), allocatable :: wf(:,:)
complex(4), allocatable :: wfp(:)
integer ntr(3),i,nrxyz(3),nrtot
integer i1,i2,i3,ir,n,m
integer ik
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: bcoeff(:,:,:,:,:,:,:)
complex(8), allocatable :: ccoeff(:,:,:)
complex(8), allocatable :: wfsvitloc(:,:,:,:)
complex(8), allocatable :: wfsvmtloc(:,:,:,:,:,:)
complex(8), allocatable :: zt2(:,:,:,:,:)
complex(8), allocatable :: zt3(:,:,:,:)
character*40 fname
real(8) x(2),alph
logical, parameter :: wf3d=.true.
logical, parameter :: wfprod=.true.
integer tlim(2,3)
integer nwfplot,firstwf
integer, external :: ikglob

call init0
call init1

! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! generate the core wavefunctions and densities
call gencore
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr

call geturf

zero3d(:)=(/8.00990d0,   4.89892d0,   7.96820d0/)

! Cartesian coordinates of boundary and origin
bound3d(:,1)=(/12.d0,0.d0,0.d0/)
bound3d(:,2)=(/0.d0,12.d0,0.d0/)
bound3d(:,3)=(/0.d0,0.d0,12.d0/)
orig3d(:)=zero3d(:)-(bound3d(:,1)+bound3d(:,2)+bound3d(:,3))/2.d0

bound2d(:,1)=(/20.d0,0.d0,0.d0/)
bound2d(:,2)=(/0.d0,20.d0,0.d0/)
orig2d(:)=(/-8.d0,-8.d0, 0.d0/)

nrxyz(:)=(/120,120,120/)

nwfplot=5
firstwf=6

if (wf3d) then
  nrtot=nrxyz(1)*nrxyz(2)*nrxyz(3)
else
  nrtot=nrxyz(1)*nrxyz(2)
endif
if (iproc.eq.0) then
  allocate(wf(nwfplot,nrtot))
  allocate(wfp(nrtot))
endif

! find the translation limits
tlim(1,:)=1000
tlim(2,:)=-1000
if (wf3d) then
  do i1=0,1
    do i2=0,1
      do i3=0,1
        r(:)=orig3d(:)+i1*bound3d(:,1)+i2*bound3d(:,2)+i3*bound3d(:,3)
        call getntr(r,ntr,vr0l)
        do i=1,3
          tlim(1,i)=min(ntr(i),tlim(1,i))
          tlim(2,i)=max(ntr(i),tlim(2,i))
        enddo
      enddo
    enddo
  enddo
else
  do i1=0,1
    do i2=0,1
      r(:)=orig2d(:)+i1*bound2d(:,1)+i2*bound2d(:,2)
      call getntr(r,ntr,vr0l)
      do i=1,3
        tlim(1,i)=min(ntr(i),tlim(1,i))
        tlim(2,i)=max(ntr(i),tlim(2,i))
      enddo
    enddo
  enddo
endif  
tlim(1,:)=tlim(1,:)-1
tlim(2,:)=tlim(2,:)+1
if (iproc.eq.0) then
  write(*,'("Translation limits:",3(4x,2I4))')tlim(:,1),tlim(:,2),tlim(:,3)
endif

allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wfsvmtloc(lmmaxvr,nrfmax,natmtot,nstsv,nspinor,nkptloc(iproc)))
allocate(wfsvitloc(ngkmax,nstsv,nspinor,nkptloc(iproc)))
! read and transform eigen-vectors
do ik=1,nkptloc(0)
  do i=0,nproc-1
    if (iproc.eq.i.and.ik.le.nkptloc(iproc)) then
      call getevecfv(vkl(1,ikglob(ik)),vgkl(1,1,1,ik),evecfv)
      call getevecsv(vkl(1,ikglob(ik)),evecsv)
      call getwfc(ikglob(ik),wann_c(1,1,1,ik))
    endif
    call barrier(comm_world)
  enddo !i
  if (ik.le.nkptloc(iproc)) then
    call match(ngk(1,ikglob(ik)),gkc(1,1,ik),tpgkc(1,1,1,ik),sfacgk(1,1,1,ik),apwalm)
    call genwfsvmt(lmaxvr,lmmaxvr,ngk(1,ikglob(ik)),evecfv,evecsv,apwalm, &
      wfsvmtloc(1,1,1,1,1,ik))
    call genwfsvit(ngk(1,ikglob(ik)),evecfv,evecsv,wfsvitloc(1,1,1,ik))
  endif
enddo

allocate(zt2(lmmaxvr,nrfmax,natmtot,nwfplot,nkptloc(iproc)))
zt2=dcmplx(0.d0,0.d0)
do n=1,nwfplot
  do ik=1,nkptloc(iproc)
    do i=1,nstfv
      zt2(:,:,:,n,ik)=zt2(:,:,:,n,ik)+wann_c(n+firstwf-1,i,1,ik)*wfsvmtloc(:,:,:,i,1,ik)
    enddo
  enddo
enddo

if (iproc.eq.0) then
  allocate(bcoeff(lmmaxvr,nrfmax,natmtot,nwfplot,tlim(1,1):tlim(2,1),tlim(1,2):tlim(2,2),tlim(1,3):tlim(2,3)))
endif

allocate(zt3(lmmaxvr,nrfmax,natmtot,nwfplot))
do i1=tlim(1,1),tlim(2,1)
do i2=tlim(1,2),tlim(2,2)
do i3=tlim(1,3),tlim(2,3)
  zt3=dcmplx(0.d0,0.d0)
  t(:)=i1*avec(:,1)+i2*avec(:,2)+i3*avec(:,3)
  do n=1,nwfplot
    do ik=1,nkptloc(iproc)
      zt3(:,:,:,n)=zt3(:,:,:,n)+exp(dcmplx(0.d0,dot_product(t,vkc(:,ikglob(ik))))) * &
        zt2(:,:,:,n,ik)/nkpt
    enddo
  enddo
  call zsync(zt3,lmmaxvr*nrfmax*natmtot*nwfplot,.true.,.false.)
  if (iproc.eq.0) then
    bcoeff(:,:,:,:,i1,i2,i3)=zt3(:,:,:,:)
  endif
enddo
enddo
enddo
deallocate(zt2)
deallocate(zt3)

allocate(ccoeff(ngkmax,nwfplot,nkptloc(iproc)))
ccoeff=dcmplx(0.d0,0.d0)
do n=1,nwfplot
  do ik=1,nkptloc(iproc)
    do i=1,nstfv
      ccoeff(:,n,ik)=ccoeff(:,n,ik)+wann_c(n+firstwf-1,i,1,ik)*wfsvitloc(:,i,1,ik)/nkpt
    enddo
  enddo
enddo

if (wf3d) then
  ir=0
  do i1=0,nrxyz(1)-1
    if (iproc.eq.0) then
      write(*,'("step ",I4," out of ",I4)')i1+1,nrxyz(1)
    endif
    do i2=0,nrxyz(2)-1
      do i3=0,nrxyz(3)-1
! arbitrary r-point
        r(:)=orig3d(:)+i1*bound3d(:,1)/nrxyz(1)+ &
                     i2*bound3d(:,2)/nrxyz(2)+ &
                     i3*bound3d(:,3)/nrxyz(3)
        ir=ir+1
        call wann_val(r,wf(:,ir),nwfplot,tlim,bcoeff,ccoeff)
      enddo
    enddo
  enddo
else
  ir=0
  do i1=0,nrxyz(1)-1
    write(*,'("step ",I4," out of ",I4)')i1+1,nrxyz(1)
    do i2=0,nrxyz(2)-1
! arbitrary r-point
      r(:)=orig2d(:)+i1*bound2d(:,1)/nrxyz(1)+ &
                     i2*bound2d(:,2)/nrxyz(2)
      ir=ir+1
      call wann_val(r,wf(:,ir),nwfplot,tlim,bcoeff,ccoeff)
    enddo
  enddo
endif

if (iproc.eq.0) then

if (.not.wf3d) then
  x(1)=sqrt(bound2d(1,1)**2+bound2d(2,1)**2+bound2d(3,1)**2)
  x(2)=sqrt(bound2d(1,2)**2+bound2d(2,2)**2+bound2d(3,2)**2)
  alph=dot_product(bound2d(:,1),bound2d(:,2))/x(1)/x(2)
  alph=acos(alph)
endif

do n=1,nwfplot
  write(fname,'("wf_",I2.2,".dx")')n
  open(70,file=trim(fname),status='replace',form='formatted')
  if (wf3d) then
    write(70,400)nrxyz(1),nrxyz(2),nrxyz(3)
    write(70,402)orig3d(:)
    do i=1,3
      write(70,404)bound3d(:,i)/nrxyz(i)
    enddo
    write(70,406)nrxyz(1),nrxyz(2),nrxyz(3)
  else
    write(70,500)nrxyz(1),nrxyz(2)
    write(70,502)(/0.d0, 0.d0/)
    write(70,504)(/x(1),0.d0/)/nrxyz(1)
    write(70,504)(/x(2)*cos(alph),x(2)*sin(alph)/)/nrxyz(2)
    write(70,506)nrxyz(1),nrxyz(2)
  endif
  write(70,408)2,nrtot
  write(70,410)(abs(wf(n,ir)),abs(wf(n,ir))**2,ir=1,nrtot)
  write(70,412)
  close(70)
enddo

if (wfprod) then
  do n=1,nwfplot
    do m=n,nwfplot
      write(fname,'("wf_prod_",I2.2,"x",I2.2,".dx")')n,m
      do ir=1,nrtot
        wfp(ir)=wf(n,ir)*conjg(wf(m,ir))
      enddo
      open(70,file=trim(fname),status='replace',form='formatted')
      if (wf3d) then
        write(70,400)nrxyz(1),nrxyz(2),nrxyz(3)
        write(70,402)orig3d(:)
        do i=1,3
          write(70,404)bound3d(:,i)/nrxyz(i)
        enddo
        write(70,406)nrxyz(1),nrxyz(2),nrxyz(3)
      else
        write(70,500)nrxyz(1),nrxyz(2)
        write(70,502)(/0.d0, 0.d0/)
        write(70,504)(/x(1),0.d0/)/nrxyz(1)
        write(70,504)(/x(2)*cos(alph),x(2)*sin(alph)/)/nrxyz(2)
        write(70,506)nrxyz(1),nrxyz(2)
      endif
      write(70,408)1,nrtot
      write(70,410)(abs(wfp(ir)),ir=1,nrtot)
      write(70,412)
      close(70)
    enddo
  enddo
endif

endif
      
 400  format('object 1 class gridpositions counts ', 3i4)
 402  format('origin ', 3f12.6)
 404  format('delta ', 3f12.6)
 406  format('object 2 class gridconnections counts ', 3i4)
 408  format('object 3 class array type float rank 1 shape',i3, &
              ' items ',i8,' data follows')
 410  format(f10.6)
 412  format('object "regular positions regular connections" class field', &
            /'component "positions" value 1',   &
	    /'component "connections" value 2', &
	    /'component "data" value 3',        &
	    /'end')

 500  format('object 1 class gridpositions counts ', 2i4)
 502  format('origin ', 2f12.6)
 504  format('delta ', 2f12.6)
 506  format('object 2 class gridconnections counts ', 2i4)







return
end


subroutine wann_val(r,val,nwfplot,tlim,bcoeff,ccoeff)
use modmain
implicit none
! arguments
real(8), intent(in) :: r(3)
integer, intent(in) :: nwfplot
integer, intent(in) :: tlim(2,3)
complex(8), intent(in) :: bcoeff(lmmaxvr,nrfmax,natmtot,nwfplot,tlim(1,1):tlim(2,1),tlim(1,2):tlim(2,2),tlim(1,3):tlim(2,3))
complex(8), intent(in) :: ccoeff(ngkmax,nwfplot,nkptloc(iproc))
complex(4), intent(out) :: val(nwfplot)

integer ivec,ntr(3),n
real(8) t(3),d,r1(3),pos(3),v1(3),t1,rmt2
integer i1,i2,i3,np2,is,ia,ias,ir,ir0,io,l,m,lm,i,j,ig,ik
real(8) ya(nprad),c(nprad),tp(2),t2
complex(8) ylm(lmmaxvr)
complex(8) zt1(nwfplot),zt2(nwfplot)
logical l1
real(8), external :: polynom 
integer, external :: ikglob

! reduce r-point to primitive cell and find corresponding translation vector
do ivec=1,3
  d=dot_product(r,bvec(:,ivec))/twopi
  ntr(ivec)=floor(d)
enddo

t(:)=ntr(1)*avec(:,1)+ntr(2)*avec(:,2)+ntr(3)*avec(:,3)

r1(:)=r(:)-t(:)

np2=nprad/2
zt1=dcmplx(0.d0,0.d0)
zt2=dcmplx(0.d0,0.d0)

l1=.false.
if (iproc.eq.0) then
! check if point is in a muffin-tin
do is=1,nspecies
  rmt2=rmt(is)**2
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do i1=-1,1
    do i2=-1,1
    do i3=-1,1
      pos(:)=atposc(:,ia,is)+i1*avec(:,1)+i2*avec(:,2)+i3*avec(:,3)
      v1(:)=r1(:)-pos(:)
      t1=v1(1)**2+v1(2)**2+v1(3)**2
      if (t1.lt.rmt2) then
        call sphcrd(v1,t1,tp)
        call genylm(lmaxvr,tp,ylm)
        do ir=1,nrmt(is)
          if (spr(ir,is).ge.t1) then
            if (ir.le.np2) then
              ir0=1
            else if (ir.gt.nrmt(is)-np2) then
              ir0=nrmt(is)-nprad+1
            else
              ir0=ir-np2
            end if
            t1=max(t1,spr(1,is))
            do io=1,nrfmax
              do l=0,lmaxvr
                do j=1,nprad
                  i=ir0+j-1
                  ya(j)=urf(i,l,io,ias)
                end do
                t2=polynom(0,nprad,spr(ir0,is),ya,c,t1)
                do m=-l,l
                  lm=idxlm(l,m)
                  do n=1,nwfplot
                    zt1(n)=zt1(n)+t2*ylm(lm)*bcoeff(lm,io,ias,n,ntr(1)+i1,ntr(2)+i2,ntr(3)+i3)
                  enddo
                end do !m
              end do  !l
            enddo !io
            l1=.true.
            goto 10
          end if
        end do !ir
      end if
    end do
    end do
    end do
  end do
end do
10 continue
endif !iproc.eq.0
call lsync(l1,1,.true.)
! otherwise use interstitial function
if (.not.l1) then
  do ik=1,nkptloc(iproc)
    do ig=1,ngk(1,ikglob(ik))
      t1=vgkc(1,ig,1,ik)*r(1)+vgkc(2,ig,1,ik)*r(2)+vgkc(3,ig,1,ik)*r(3)
      do n=1,nwfplot
        zt2(n)=zt2(n)+cmplx(cos(t1),sin(t1),8)*ccoeff(ig,n,ik)/sqrt(omega)
      enddo
    end do
  end do
  call zsync(zt2,nwfplot,.true.,.false.)
endif
if (iproc.eq.0) then
  val(:)=zt1(:)+zt2(:)
endif
return
end


!subroutine getntr(v,ntr)
!use modmain
!implicit none
!real(8), intent(in) :: v(3)
!integer, intent(out) :: ntr(3)
!real(8) t
!integer ivec
!
!do ivec=1,3
!  t=dot_product(v,bvec(:,ivec))/twopi
!  ntr(ivec)=floor(t)
!enddo
!
!return
!end

subroutine getntr(vrc,ntr,vr0l)
use modmain
implicit none
real(8), intent(in) :: vrc(3)
integer, intent(out) :: ntr(3)
real(8), intent(out) :: vr0l(3)
real(8) a(3,3)
real(8) b(3)
integer i,j,ipiv(3)
real(8) work(24)
integer lwork

do j=1,3
  do i=1,j
    a(i,j)=dot_product(avec(:,i),avec(:,j))
  enddo
  b(j)=dot_product(avec(:,j),vrc(:))
enddo
lwork=24
call dsysv('U',3,1,a,3,ipiv,b,3,work,lwork,i)
do i=1,3
  ntr(i)=floor(b(i))
  vr0l(i)=b(i)-ntr(i)
enddo

return
end

subroutine putwfc(ik,wfcp)
use modmain
implicit none
integer, intent(in) :: ik
complex(8), intent(in) :: wfcp(wann_nmax,nstfv,wann_nspin)
integer recl

inquire(iolength=recl)ik,wann_nmax,nstfv,wann_nspin,wfcp
open(70,file='WANN_C.OUT',action='WRITE',form='UNFORMATTED', &
  access='DIRECT',recl=recl)
write(70,rec=ik)ik,wann_nmax,nstfv,wann_nspin,wfcp
close(70)

return
end

subroutine getwfc(ik,wfcp)
use modmain
implicit none
integer, intent(in) :: ik
complex(8), intent(out) :: wfcp(wann_nmax,nstfv,wann_nspin)
integer ik_,wann_nmax_,nstfv_,wann_nspin_,recl

inquire(iolength=recl)ik_,wann_nmax_,nstfv_,wann_nspin_,wfcp
open(70,file='WANN_C.OUT',action='READ',form='UNFORMATTED', &
  access='DIRECT',recl=recl)
read(70,rec=ik)ik_,wann_nmax_,nstfv_,wann_nspin_,wfcp
close(70)
if (ik_.ne.ik.or.wann_nmax_.ne.wann_nmax.or.nstfv_.ne.nstfv.or. &
  wann_nspin_.ne.wann_nspin) then
  write(*,*)
  write(*,'("Error(getwfc): wrong dimensions")')
  write(*,*)
  call pstop
endif

return
end

subroutine putwann(ik)
use modmain
implicit none
integer, intent(in) :: ik
integer recl
integer, external :: ikglob

inquire(iolength=recl)ik,wann_nmax,wann_nspin,lmmaxvr,nrfmax,natmtot,ngkmax, &
  wann_unkmt(:,:,:,:,:,ik),wann_unkit(:,:,:,ik)
open(70,file='WANN_UNK.OUT',action='WRITE',form='UNFORMATTED', &
  access='DIRECT',recl=recl)
write(70,rec=ikglob(ik))ikglob(ik),wann_nmax,wann_nspin,lmmaxvr,nrfmax,natmtot,ngkmax, &
  wann_unkmt(:,:,:,:,:,ik),wann_unkit(:,:,:,ik)
close(70)

return
end

subroutine getwann(ik)
use modmain
implicit none
integer, intent(in) :: ik
integer ik_,wann_nmax_,wann_nspin_,lmmaxvr_,nrfmax_,natmtot_,ngkmax_,recl
integer, external :: ikglob

inquire(iolength=recl)ik_,wann_nmax_,wann_nspin_,lmmaxvr_,nrfmax_,natmtot_,ngkmax_, &
  wann_unkmt(:,:,:,:,:,ik),wann_unkit(:,:,:,ik)
open(70,file='WANN_UNK.OUT',action='READ',form='UNFORMATTED',access='DIRECT',recl=recl)
read(70,rec=ikglob(ik))ik_,wann_nmax_,wann_nspin_,lmmaxvr_,nrfmax_,natmtot_,ngkmax_, &
  wann_unkmt(:,:,:,:,:,ik),wann_unkit(:,:,:,ik)
close(70)
if (ik_.ne.ikglob(ik).or.wann_nmax_.ne.wann_nmax.or.wann_nspin_.ne.wann_nspin.or. &
  lmmaxvr_.ne.lmmaxvr.or.nrfmax_.ne.nrfmax.or.natmtot_.ne.natmtot.or.ngkmax_.ne.ngkmax) then
  write(*,*)
  write(*,'("Error(getwann): wrong dimensions")')
  write(*,'("  proc : ",I4)')iproc
  write(*,'("  ik_ : ",I4," ik : ",I4)')ik_,ikglob(ik)
  write(*,'("  wann_nmax_ : ",I4," wann_nmax : ",I4)')wann_nmax_,wann_nmax
  write(*,'("  wann_nspin_ : ",I4," wann_nspin : ",I4)')wann_nspin_,wann_nspin
  write(*,'("  lmmaxvr_ : ",I4," lmmaxvr : ",I4)')lmmaxvr_,lmmaxvr
  write(*,'("  nrfmax_ : ",I4," nrfmax : ",I4)')nrfmax_,nrfmax
  write(*,'("  natmtot_ : ",I4," natmtot : ",I4)')natmtot_,natmtot
  write(*,'("  ngkmax_ : ",I4," ngkmax : ",I4)')ngkmax_,ngkmax
  write(*,*)
  call pstop
endif

return
end



