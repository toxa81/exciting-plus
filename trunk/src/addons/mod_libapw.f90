module mod_libapw

integer nrfmtmax
integer ordrfmtmax
integer, allocatable :: ordrfmt(:,:)
integer, allocatable :: nrfmt(:)
real(8), allocatable :: hmltrad(:,:,:,:)
real(8), allocatable :: beffrad(:,:,:,:,:)
real(8), allocatable :: ovlprad(:,:,:,:)
real(8), allocatable :: beffir(:,:)

contains

subroutine libapw_init 
use modmain
implicit none
integer is,l,io,ilo,ik,ikloc

#ifdef _LIBAPW_
if (allocated(nrfmt)) deallocate(nrfmt)
allocate(nrfmt(nspecies))
if (allocated(ordrfmt)) deallocate(ordrfmt)
allocate(ordrfmt(0:lmaxapw,nspecies))
nrfmt=0
ordrfmt=0
do is=1,nspecies
  do l=0,lmaxapw
    do io=1,apword(l,is)
      nrfmt(is)=nrfmt(is)+1
      ordrfmt(l,is)=ordrfmt(l,is)+1
    enddo
  enddo
  do ilo=1,nlorb(is)
    nrfmt(is)=nrfmt(is)+1
    ordrfmt(lorbl(ilo,is),is)=ordrfmt(lorbl(ilo,is),is)+1
  enddo
enddo
nrfmtmax=maxval(nrfmt)
ordrfmtmax=maxval(ordrfmt)
if (allocated(hmltrad)) deallocate(hmltrad)
allocate(hmltrad(lmmaxvr,nrfmtmax,nrfmtmax,natmtot))
hmltrad=0.d0
if (allocated(beffrad)) deallocate(beffrad)
allocate(beffrad(lmmaxvr,nrfmtmax,nrfmtmax,natmtot,ndmag))
hmltrad=0.d0
if (allocated(ovlprad)) deallocate(ovlprad)
allocate(ovlprad(0:lmaxapw,ordrfmtmax,ordrfmtmax,natmtot))
ovlprad=0.d0
if (allocated(beffir)) deallocate(beffir)
allocate(beffir(ngrtot,ndmag))
beffir=0.d0
call lapw_load_global(natmtot,nspecies,lmaxvr,lmaxapw,apwordmax,nrmtmax,&
  &ngkmax,ngvec,ngrtot,nlomax,ias2is,intgv,ivg,ivgig,ngrid,igfft,cfunir, &
  &cfunig,gntyry,nstfv,nstsv,nmatmax,nrfmtmax,ordrfmtmax,evaltol,spinpol,ndmag)
do is=1,nspecies
  call lapw_load_species(is,nlorb(is),lorbl(1,is),apword(0,is),rmt(is),nrmt(is))
enddo
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
  call lapw_load_kpoint(ngk(1,ik),igkig(1,1,ikloc),vgkc(1,1,1,ikloc))
enddo
call lapw_init
#endif
return
end subroutine

#ifdef _LIBAPW_
subroutine libapw_seceqn_init
use modmain
implicit none
integer ir,is,ia,ias,l1,l2,io1,io2,ilo1,ilo2,i1,i2,nr,i
integer l1tmp(0:lmaxapw),l2tmp,lm
real(8) cb,t1
real(8), allocatable :: rfmt(:,:)
real(8), allocatable :: bmt(:,:,:)
real(8) r2(nrmtmax),fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)
!
allocate(rfmt(nrmtmax,nrfmtmax))
allocate(bmt(lmmaxvr,nrmtmax,ndmag))
cb=gfacte/(4.d0*solsc)
do ias=1,natmtot
  is=ias2is(ias)
  ia=ias2ia(ias)
  nr=nrmt(is)
  do ir=1,nr
    r2(ir)=spr(ir,is)**2
  end do
  i2=0
  do l2=0,lmaxapw
    do io2=1,apword(l2,is)
      i2=i2+1
      i1=0
      do l1=0,lmaxapw
        do io1=1,apword(l1,is)
          i1=i1+1
          hmltrad(:,i1,i2,ias)=haa(:,io1,l1,io2,l2,ias)
        enddo
      enddo
    enddo !io2
  enddo !l2
  do ilo2=1,nlorb(is)
    i2=i2+1
    i1=0
    do l1=0,lmaxapw
      do io1=1,apword(l1,is)
        i1=i1+1
        hmltrad(:,i1,i2,ias)=hloa(:,ilo2,io1,l1,ias)
      enddo
    enddo
    do ilo1=1,nlorb(is)
      i1=i1+1
      hmltrad(:,i1,i2,ias)=hlolo(:,ilo1,ilo2,ias)
    enddo
  enddo !ilo2
! overlap integrals
  do l1=0,lmaxapw
    do io1=1,apword(l1,is)
      ovlprad(l1,io1,io1,ias)=1.d0
    enddo
  enddo
  l1tmp=0
  do ilo1=1,nlorb(is)
    l1=lorbl(ilo1,is)
    l1tmp(l1)=l1tmp(l1)+1
    do io1=1,apword(l1,is)
      ovlprad(l1,io1,apword(l1,is)+l1tmp(l1),ias)=oalo(io1,ilo1,ias)
      ovlprad(l1,apword(l1,is)+l1tmp(l1),io1,ias)=oalo(io1,ilo1,ias)
    enddo
    l2tmp=0
    do ilo2=1,nlorb(is)
      if (l1.eq.lorbl(ilo2,is)) then
        l2tmp=l2tmp+1
        ovlprad(l1,apword(l1,is)+l1tmp(l1),apword(l1,is)+l2tmp,ias)=ololo(ilo1,ilo2,ias)
        ovlprad(l1,apword(l1,is)+l2tmp,apword(l1,is)+l1tmp(l1),ias)=ololo(ilo2,ilo1,ias)
      endif
    enddo
  enddo
! collect radial functions
  i1=0
  do l1=0,lmaxapw
    do io1=1,apword(l1,is)
      i1=i1+1
      rfmt(:,i1)=apwfr(:,1,io1,l1,ias)
    enddo
  enddo
  do ilo1=1,nlorb(is)
    i1=i1+1
    rfmt(:,i1)=lofr(:,1,ilo1,ias)
  enddo
! generate radial integals for magnetic field
  if (spinpol) then
    bmt=0.d0
    ! z
    bmt(:,:,1)=bxcmt(:,:,ias,ndmag)
    t1=cb*(bfcmt(3,ia,is)+bfieldc(3))
    bmt(1,:,1)=bmt(1,:,1)+t1/y00
    if (ndmag.eq.3) then
      ! x
      bmt(:,:,2)=bxcmt(:,:,ias,1)
      t1=cb*(bfcmt(1,ia,is)+bfieldc(1))
      bmt(1,:,2)=bmt(1,:,2)+t1/y00
      ! y
      bmt(:,:,3)=bxcmt(:,:,ias,2)
      t1=cb*(bfcmt(2,ia,is)+bfieldc(2))
      bmt(1,:,3)=bmt(1,:,3)+t1/y00
    endif
    do i=1,ndmag
      do i1=1,nrfmt(is)
        do i2=1,nrfmt(is)
          do lm=1,lmmaxvr
            do ir=1,nr
              fr(ir)=rfmt(ir,i1)*rfmt(ir,i2)*r2(ir)*bmt(lm,ir,i)
            enddo
            call fderiv(-1,nr,spr(:,is),fr,gr,cf)
            beffrad(lm,i1,i2,ias,i)=gr(nr)
          enddo !lm
        enddo
      enddo
    enddo !i
  endif ! spinpol
enddo !ias
deallocate(rfmt,bmt)
if (spinpol) then
  ! z
  do ir=1,ngrtot
    beffir(ir,1)=bxcir(ir,ndmag)+cb*bfieldc(3)
  enddo
  if (ndmag.eq.3) then
    ! x
    do ir=1,ngrtot
      beffir(ir,2)=bxcir(ir,1)+cb*bfieldc(1)
    enddo
    ! y
    do ir=1,ngrtot
      beffir(ir,3)=bxcir(ir,2)+cb*bfieldc(2)
    enddo
  endif
endif
call lapw_seceqn_init(hmltrad,ovlprad,beffrad,apwfr,apwdfr,beffir,veffig) 
return
end subroutine
#endif

end module
