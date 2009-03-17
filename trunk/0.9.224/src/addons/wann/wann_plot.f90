subroutine wann_plot
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none

real(8) orig(3)
complex(4), allocatable :: wf(:,:)
complex(4), allocatable :: wfp(:)
integer i,nrtot
integer i1,i2,i3,ir,n,m
real(8), allocatable :: vr(:,:)
real(8), allocatable :: veff(:)
complex(8), allocatable :: zfft_vir(:)
character*40 fname
real(8) x(2),alph
logical, parameter :: wfprod=.false.
integer ikloc
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

do i=0,nproc-1
  if (iproc.eq.i) then
    do ikloc=1,nkptloc(iproc)
      call getwann(ikloc)
    end do
  end if
  call barrier(comm_world)
end do

if (task.eq.361) then
  nrtot=nrxyz(1)
  nrxyz(2)=1
  nrxyz(3)=1
  orig(:)=zero3d(:)-(bound3d(:,1))/2.d0
endif
if (task.eq.362) then
  nrtot=nrxyz(1)*nrxyz(2)
  nrxyz(3)=1
  orig(:)=zero3d(:)-(bound3d(:,1)+bound3d(:,2))/2.d0
endif
if (task.eq.363) then
  nrtot=nrxyz(1)*nrxyz(2)*nrxyz(3)
  orig(:)=zero3d(:)-(bound3d(:,1)+bound3d(:,2)+bound3d(:,3))/2.d0
endif

if (iproc.eq.0) then
  allocate(wf(nwfplot,nrtot))
  allocate(wfp(nrtot))
  allocate(veff(nrtot))
endif
allocate(vr(3,nrtot))

! make (1,2,3)D-grid of r-points
ir=0
do i1=0,nrxyz(1)-1
  do i2=0,nrxyz(2)-1
    do i3=0,nrxyz(3)-1
      ir=ir+1
      vr(:,ir)=orig(:)+i1*bound3d(:,1)/nrxyz(1)+ &
                       i2*bound3d(:,2)/nrxyz(2)+ &
                       i3*bound3d(:,3)/nrxyz(3)
    enddo
  enddo
enddo
if (ir.ne.nrtot) then
  write(*,*)'Error in r-mesh'
  call pstop
endif

! Fourier transform potential to G-space
allocate(zfft_vir(ngrtot))
zfft_vir(:)=veffir(:)
call zfftifc(3,ngrid,-1,zfft_vir)

do ir=1,nrtot
  if (mod(ir,nrxyz(2)*nrxyz(3)).eq.0.and.iproc.eq.0) then
    write(*,*)'r-point : ',ir,' out of ',nrtot
  endif
  call wann_val(vr(1,ir),wf(1,ir))
  if (iwfv.ne.0) call f_veff_p(vr(1,ir),zfft_vir,veff(ir))
enddo

if (iproc.eq.0) then
  write(*,*)'MT part:',timer(1,2)
  write(*,*)'IT part:',timer(2,2)

  if (task.eq.362) then
    x(1)=sqrt(bound3d(1,1)**2+bound3d(2,1)**2+bound3d(3,1)**2)
    x(2)=sqrt(bound3d(1,2)**2+bound3d(2,2)**2+bound3d(3,2)**2)
    alph=dot_product(bound3d(:,1),bound3d(:,2))/x(1)/x(2)
    alph=acos(alph)
  endif
  
  do n=1,nwfplot
    write(fname,'("wf_",I3.3,".dx")')n+firstwf-1
    open(70,file=trim(fname),status='replace',form='formatted')
    if (task.eq.363) then
      write(70,400)nrxyz(1),nrxyz(2),nrxyz(3)
      write(70,402)orig(:)
      do i=1,3
        write(70,404)bound3d(:,i)/nrxyz(i)
      enddo
      write(70,406)nrxyz(1),nrxyz(2),nrxyz(3)
    endif
    if (task.eq.362) then
      write(70,500)nrxyz(1),nrxyz(2)
      write(70,502)(/0.d0, 0.d0/)
      write(70,504)(/x(1),0.d0/)/nrxyz(1)
      write(70,504)(/x(2)*cos(alph),x(2)*sin(alph)/)/nrxyz(2)
      write(70,506)nrxyz(1),nrxyz(2)
    endif
    write(70,408)2,nrtot
    write(70,'(4G18.10)')(abs(wf(n,ir)),abs(wf(n,ir))**2,ir=1,nrtot)
    write(70,412)
    close(70)
  enddo
  
  write(fname,'("veff.dx")')
  open(70,file=trim(fname),status='replace',form='formatted')
  if (task.eq.363) then
    write(70,400)nrxyz(1),nrxyz(2),nrxyz(3)
    write(70,402)orig(:)
    do i=1,3
      write(70,404)bound3d(:,i)/nrxyz(i)
    enddo
    write(70,406)nrxyz(1),nrxyz(2),nrxyz(3)
  endif
  if (task.eq.362) then
    write(70,500)nrxyz(1),nrxyz(2)
    write(70,502)(/0.d0, 0.d0/)
    write(70,504)(/x(1),0.d0/)/nrxyz(1)
    write(70,504)(/x(2)*cos(alph),x(2)*sin(alph)/)/nrxyz(2)
    write(70,506)nrxyz(1),nrxyz(2)
  endif
  write(70,408)1,nrtot
  write(70,'(4G18.10)')(veff(ir),ir=1,nrtot)
  write(70,412)
  close(70)
  
!  if (wfprod) then
!    do n=1,nwfplot
!      do m=n,nwfplot
!        write(fname,'("wf_prod_",I2.2,"x",I2.2,".dx")')n,m
!        do ir=1,nrtot
!          wfp(ir)=wf(n,ir)*conjg(wf(m,ir))
!        enddo
!        open(70,file=trim(fname),status='replace',form='formatted')
!        if (wf3d) then
!          write(70,400)nrxyz(1),nrxyz(2),nrxyz(3)
!          write(70,402)orig3d(:)
!          do i=1,3
!            write(70,404)bound3d(:,i)/nrxyz(i)
!          enddo
!          write(70,406)nrxyz(1),nrxyz(2),nrxyz(3)
!        else
!          write(70,500)nrxyz(1),nrxyz(2)
!          write(70,502)(/0.d0, 0.d0/)
!          write(70,504)(/x(1),0.d0/)/nrxyz(1)
!          write(70,504)(/x(2)*cos(alph),x(2)*sin(alph)/)/nrxyz(2)
!          write(70,506)nrxyz(1),nrxyz(2)
!        endif
!        write(70,408)1,nrtot
!        write(70,410)(abs(wfp(ir)),ir=1,nrtot)
!        write(70,412)
!        close(70)
!      enddo
!    enddo
!  endif

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


subroutine wann_val(r,val)
use modmain
implicit none
! arguments
real(8), intent(in) :: r(3)
complex(4), intent(out) :: val(nwfplot)

integer ntr(3),n
integer is,ia,ias,ir0,io,lm,ig,i,j,l
real(8) ya(nprad),c(nprad),tp(2)
real(8) vr0(3),t1,r0,tr(3)
real(8) ur(0:lmaxvr,nrfmax)
complex(8) ylm(lmmaxvr)
complex(8) zt1(nwfplot),zt2(nwfplot),zt3
integer ikloc
real(8), external :: polynom 
integer, external :: ikglob
logical, external :: vrinmt

zt1=dcmplx(0.d0,0.d0)
zt2=dcmplx(0.d0,0.d0)
if (vrinmt(r,is,ia,ntr,vr0,ir0,r0)) then
  call timer_start(1)
  ias=idxas(ia,is)
  call sphcrd(vr0,t1,tp)
  call genylm(lmaxvr,tp,ylm)
  tr(:)=ntr(1)*avec(:,1)+ntr(2)*avec(:,2)+ntr(3)*avec(:,3)
  do io=1,nrfmax
    do l=0,lmaxvr
      do j=1,nprad
        i=ir0+j-1
        ya(j)=urf(i,l,io,ias)
      end do
      ur(l,io)=polynom(0,nprad,spr(ir0,is),ya,c,r0)
    enddo !l
  enddo !io
  do ikloc=1,nkptloc(iproc)
    zt3=exp(zi*dot_product(vkc(:,ikglob(ikloc)),tr(:)))
    do n=1,nwfplot
      do io=1,nrfmax
        do lm=1,lmmaxvr
          zt1(n)=zt1(n)+zt3*wann_unkmt(lm,io,ias,n+firstwf-1,1,ikloc)* &
            ur(lm2l(lm),io)*ylm(lm)
        enddo
      enddo
    enddo
  enddo
  call zsync(zt1,nwfplot,.true.,.false.)
  call timer_stop(1)
else
  call timer_start(2)
  do ikloc=1,nkptloc(iproc)
    do ig=1,ngk(1,ikglob(ikloc))
      zt3=exp(zi*dot_product(vgkc(:,ig,1,ikloc),r(:)))/sqrt(omega)
      do n=1,nwfplot
        zt2(n)=zt2(n)+zt3*wann_unkit(ig,n+firstwf-1,1,ikloc)
      enddo
    end do
  end do
  call zsync(zt2,nwfplot,.true.,.false.)
  call timer_stop(2)
endif
if (iproc.eq.0) then
  val(:)=(zt1(:)+zt2(:))/nkpt
endif
call barrier(comm_world)
return
end



logical function vrinmt(vr,is,ia,tr,vr0,ir0,r0)
use modmain
implicit none
! arguments
real(8), intent(in) :: vr(3)
integer, intent(out) :: is
integer, intent(out) :: ia
integer, intent(out) :: tr(3)
real(8), intent(out) :: vr0(3)
integer, intent(out) :: ir0
real(8), intent(out) :: r0
! local variables
real(8) vr0l(3),r1(3),rmt2,pos(3)
integer i1,i2,i3,ir,np2

vrinmt=.false.
call getntr(vr,tr,vr0l)
r1(:)=vr0l(1)*avec(:,1)+vr0l(2)*avec(:,2)+vr0l(3)*avec(:,3)

np2=nprad/2
do is=1,nspecies
  rmt2=rmt(is)**2
  do ia=1,natoms(is)
    do i1=-1,1
    do i2=-1,1
    do i3=-1,1
      pos(:)=atposc(:,ia,is)+i1*avec(:,1)+i2*avec(:,2)+i3*avec(:,3)
      vr0(:)=r1(:)-pos(:)
      r0=vr0(1)**2+vr0(2)**2+vr0(3)**2
      if (r0.lt.rmt2) then
        r0=sqrt(r0)
        do ir=1,nrmt(is)
          if (spr(ir,is).ge.r0) then
            if (ir.le.np2) then
              ir0=1
            else if (ir.gt.nrmt(is)-np2) then
              ir0=nrmt(is)-nprad+1
            else
              ir0=ir-np2
            end if
            r0=max(r0,spr(1,is))
            tr(1)=tr(1)+i1
            tr(2)=tr(2)+i2
            tr(3)=tr(3)+i3
            vrinmt=.true.
            return
          end if
        end do !ir
      end if
    end do
    end do
    end do
  end do
end do
return
end

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

! r=r0+T=r0+n1*a1+n2*a2+n3*a3=f1*a1+f2*a2+f3*a3
! system of linear equtions for f1,f2,f3
! r*a1=f1*a1*a1+f2*a2*a1+f3*a3*a1
! r*a2=f1*a1*a2+f2*a2*a2+f3*a3*a2
! r*a3=f1*a1*a3+f2*a2*a3+f3*a3*a3
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



