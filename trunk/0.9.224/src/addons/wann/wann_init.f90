subroutine wann_init
use modmain
implicit none

integer i,j,n,lm,ispn,iwgrp,itype

if (nspnfv.eq.2) then
  write(*,*)
  write(*,'("Error(wann_init): can''t do Wannier + spin spirals")')
  write(*,*)
  call pstop
endif

if (ncmag) then
  write(*,*)
  write(*,'("Error(wann_init): can''t do Wannier + non-collinear magnetism")')
  write(*,*)
  call pstop
endif

nwann(:)=0
do ispn=1,wann_nspin
  do i=1,wann_natom
    iwgrp=wann_iatom(1+ispn,i)
    nwann(ispn)=nwann(ispn)+wann_iorbgrp(0,1,iwgrp)
  enddo
enddo
wann_nmax=maxval(nwann)
if (allocated(iwann)) deallocate(iwann)
allocate(iwann(wann_nmax,wann_nspin,4))
if (allocated(iasiwann)) deallocate(iasiwann)
allocate(iasiwann(natmtot,lmmaxlu,nspinor))
iasiwann=-1

do ispn=1,wann_nspin
  n=0
  do i=1,wann_natom
    iwgrp=wann_iatom(1+ispn,i)
    do j=1,wann_iorbgrp(0,1,iwgrp)
      n=n+1
      lm=wann_iorbgrp(j,1,iwgrp)
      itype=wann_iorbgrp(j,2,iwgrp)
      iwann(n,ispn,1)=wann_iatom(1,i)
      iwann(n,ispn,2)=lm
      iwann(n,ispn,3)=lm2l(lm)
      iwann(n,ispn,4)=itype
      iasiwann(iwann(n,ispn,1),iwann(n,ispn,2),ispn)=n
    enddo
  enddo !i
enddo !ispn
      
open(100,file='WANNIER.OUT',form='formatted',status='replace')
write(100,'("number of WF atoms : ", I4)')wann_natom
write(100,'("number of WF spins : ", I1)')wann_nspin
write(100,'("number of WF orbital groups : ", I4)')wann_norbgrp
write(100,'("number of WF types : ", I4)')wann_ntype
write(100,*)
do ispn=1,wann_nspin
  write(100,'("spin : ",I1,"    number of WF : ",I4)')ispn,nwann(ispn)
  do n=1,nwann(ispn)
    itype=iwann(n,ispn,4)
    write(100,'("  wf : ",I4)')n
    write(100,'("    atom : ",I4)')iwann(n,ispn,1)
    write(100,'("    lm   : ",I4)')iwann(n,ispn,2)
    write(100,'("    l    : ",I4)')iwann(n,ispn,3)
    write(100,'("    type : ",I4)')iwann(n,ispn,4)
    write(100,'("    energy interval : [",F8.4,",",F8.4,"]")')wann_eint(:,itype)
    write(100,'("    band interval : from ",I4," to ",I4)')wann_nint(:,itype)
    write(100,'("    potential : ",F8.4)')wann_v(itype)
  enddo
enddo
close(100)

if (allocated(wann_c)) deallocate(wann_c)
allocate(wann_c(wann_nmax,nstfv,wann_nspin,nkptloc(iproc)))
wann_c=dcmplx(0.d0,0.d0)
if (allocated(wann_unkmt)) deallocate(wann_unkmt)
allocate(wann_unkmt(lmmaxvr,nrfmax,natmtot,wann_nmax,wann_nspin,nkptloc(iproc)))
wann_unkmt=dcmplx(0.d0,0.d0)
if (allocated(wann_unkit)) deallocate(wann_unkit)
allocate(wann_unkit(ngkmax,wann_nmax,wann_nspin,nkptloc(iproc)))
wann_unkit=dcmplx(0.d0,0.d0)

if (allocated(wann_h)) deallocate(wann_h)
allocate(wann_h(wann_nmax,wann_nmax,wann_nspin,nkpt))
wann_h=dcmplx(0.d0,0.d0)
if (allocated(wann_e)) deallocate(wann_e)
allocate(wann_e(wann_nmax,wann_nspin,nkpt))
wann_e=0.d0

if (allocated(wann_ene)) deallocate(wann_ene)
allocate(wann_ene(wann_nmax,wann_nspin))
wann_ene=0.d0
if (allocated(wann_occ)) deallocate(wann_occ)
allocate(wann_occ(wann_nmax,wann_nspin))
wann_occ=0.d0

!if (allocated(wannmt)) deallocate(wannmt)
!allocate(wannmt(lmmaxvr,nrcmtmax,natmtot,nspinor,wann_nmax,nkptloc(iproc)))
!if (allocated(wannit)) deallocate(wannit)
!allocate(wannit(ngrtot,nspinor,wann_nmax,nkptloc(iproc)))

!if (allocated(wf_p)) deallocate(wf_p)
!allocate(wf_p(3,wf_dim,wf_dim,wann_nspin,nkpt))
!wf_p=dcmplx(0.d0,0.d0)



return
end
