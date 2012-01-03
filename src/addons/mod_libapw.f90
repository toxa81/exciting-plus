module mod_libapw

integer nrfmtmax
integer, allocatable :: nrfmt(:)
real(8), allocatable :: hmltrad(:,:,:,:)
real(8), allocatable :: ovlprad(:,:,:)

contains

subroutine libapw_init 
use modmain
implicit none
integer is,l,io,ilo

#ifdef _LIBAPW_
if (allocated(nrfmt)) deallocate(nrfmt)
allocate(nrfmt(nspecies))
if (allocated(ovlprad)) deallocate(ovlprad)
allocate(ovlprad(apwordmax+nlomax,nlomax,natmtot))
ovlprad=0.d0

nrfmt=0
do is=1,nspecies
  do l=0,lmaxapw
    do io=1,apword(l,is)
      nrfmt(is)=nrfmt(is)+1
    enddo
  enddo
  do ilo=1,nlorb(is)
    nrfmt(is)=nrfmt(is)+1
  enddo
enddo
nrfmtmax=maxval(nrfmt)
if (allocated(hmltrad)) deallocate(hmltrad)
allocate(hmltrad(lmmaxvr,nrfmtmax,nrfmtmax,natmtot))
hmltrad=0.d0
call lapw_load_global(natmtot,nspecies,lmaxvr,lmaxapw,apwordmax,nrmtmax,&
  ngkmax,ngvec,ngrtot,nlomax,ias2is,intgv,ivg,ivgig,cfunig,gntyry,nstfv,&
  nstsv,nmatmax,nrfmtmax,evaltol)
do is=1,nspecies
  call lapw_load_species(is,nlorb(is),lorbl(1,is),apword(0,is),rmt(is),nrmt(is))
enddo
call lapw_init
#endif
return
end subroutine

#ifdef _LIBAPW_
subroutine libapw_radial_integrals
use modmain
implicit none
integer is,ias,l1,l2,io1,io2,ilo1,ilo2,i1,i2

do ias=1,natmtot
  is=ias2is(ias)
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
  do ilo2=1,nlomax
    ovlprad(1:apwordmax,ilo2,ias)=oalo(1:apwordmax,ilo2,ias)
    ovlprad(apwordmax+1:apwordmax+nlomax,ilo2,ias)=ololo(1:nlomax,ilo2,ias)
  enddo
enddo !ias
return
end subroutine
#endif

end module
