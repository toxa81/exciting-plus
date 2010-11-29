program main
implicit none
      
integer lmmax,natmtot,nspinor,nstfv,nstsv,nkpt,nlines
integer lm,ias,ispn,ist,ik,spin
real(4), allocatable :: bndchr(:,:,:,:,:)
real(8), allocatable :: dpp1d(:),e(:,:),w(:,:),lines(:,:)
integer, allocatable :: orb(:,:),tmp(:),wf1(:),wf2(:)
integer i,j,ist1,ist2,n,iin
real(8) scale,emin,emax,wmax,wmax1,efermi
real(8), parameter :: ha2ev = 27.21138386d0
logical wannier,l1,unitsev
integer wf_dim
real(8), allocatable :: wfc(:,:,:,:)
character*100 st
real(8) t1
integer, allocatable :: bndf(:,:)
integer nspecies,ias1,is1,is
character(256), allocatable :: spsymb(:)
integer, allocatable :: natoms(:)
integer, allocatable :: iasis(:)

      
open(50,file="BNDCHR.OUT",form="FORMATTED",status="OLD")
read(50,*)lmmax,nspecies,natmtot,nspinor,nstfv,nstsv,nkpt,nlines
read(50,*)efermi
allocate(bndchr(lmmax,natmtot,nspinor,nstsv,nkpt))
allocate(dpp1d(nkpt))
allocate(e(nstsv,nkpt))
allocate(spsymb(nspecies))
allocate(natoms(nspecies))
allocate(iasis(natmtot))
do is=1,nspecies
  read(50,*)spsymb(is)
  read(50,*)natoms(is)
enddo
do ias=1,natmtot
  read(50,*)ias1,is1
  iasis(ias1)=is1
enddo      
do ik=1,nkpt
  read(50,*)dpp1d(ik)
  read(50,*)(e(ist,ik),ist=1,nstsv)
  read(50,*)((((bndchr(lm,ias,ispn,ist,ik),lm=1,lmmax), &
                ias=1,natmtot),ispn=1,nspinor),ist=1,nstsv)
enddo
  
allocate(w(nstsv,nkpt))
allocate(orb(lmmax,natmtot))
allocate(tmp(lmmax))
      
!read(50,*)wannier
!if (wannier) then
!  read(50,*)wf_dim
!  allocate(wfc(wf_dim,nstfv,nspinor,nkpt))
!  do ik=1,nkpt
!    read(50,*)((wfc(n,i,1,ik),n=1,wf_dim),i=1,nstfv)
!  enddo
!endif
close(50)

write(*,'("Input bottom band : ")')
read(*,*)ist1
write(*,'("Input top band : ")')
read(*,*)ist2
write(*,'("Input spinor component : ")')
read(*,*)ispn

allocate(bndf(lmmax,natmtot))
bndf=1
do ik=1,nkpt
  do ias=1,natmtot
    do lm=1,lmmax
      t1=0.d0
      do ist=ist1,ist2
        t1=t1+bndchr(lm,ias,ispn,ist,ik)
      enddo
      if (t1.lt.1d-2) bndf(lm,ias)=0
    enddo
  enddo
enddo

do ias=1,natmtot
  write(*,'("Atom : ",I4)')ias
  do lm=1,lmmax
    t1=0.d0
    do ist=ist1,ist2
      do ik=1,nkpt
        t1=t1+bndchr(lm,ias,ispn,ist,ik)
      enddo    
    enddo
    write(*,'(" lm : ",I4,"  weight : ",F12.6,"  in all k : ",I1)')lm,t1,bndf(lm,ias)
  enddo
enddo
  


      
end
      
