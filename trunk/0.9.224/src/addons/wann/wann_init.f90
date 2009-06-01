subroutine wann_init
use modmain
implicit none

integer i,j,n,lm,ispn,iwgrp,itype,iatom
integer, allocatable :: iwann_tmp(:,:)

if (nspnfv.eq.2) then
  write(*,*)
  write(*,'("Error(wann_init): can''t do Wannier + spin spirals")')
  write(*,*)
  call pstop
endif

nwann=0
do i=1,wann_natom
  iwgrp=wann_iprj(2,i)
  nwann=nwann+wann_norb(iwgrp)
enddo
if (allocated(iwann)) deallocate(iwann)
allocate(iwann(4,nwann))

n=0
do i=1,wann_natom
  iatom=wann_iprj(1,i)
  iwgrp=wann_iprj(2,i)
  do j=1,wann_norb(iwgrp)
    n=n+1
    lm=wann_iorb(1,j,iwgrp)
    ispn=wann_iorb(2,j,iwgrp)
    itype=wann_iorb(3,j,iwgrp)
    iwann(1,n)=iatom
    iwann(2,n)=lm
    iwann(3,n)=ispn
    iwann(4,n)=itype
  enddo
enddo !i
! sort WF by spin (for collinear case)
if (.not.ncmag) then
  allocate(iwann_tmp(4,nwann))
  n=0
  do j=1,nwann
    if (iwann(3,j).eq.1) then
      n=n+1
      iwann_tmp(:,n)=iwann(:,j)
    endif
  enddo
  do j=1,nwann
    if (iwann(3,j).eq.2) then
      n=n+1
      iwann_tmp(:,n)=iwann(:,j)
    endif
  enddo
  iwann=iwann_tmp
  deallocate(iwann_tmp)
endif
      
open(100,file='WANNIER.OUT',form='formatted',status='replace')
write(100,'("number of WF atoms : ", I4)')wann_natom
write(100,'("number of WF orbital groups : ", I4)')wann_norbgrp
write(100,'("number of WF types : ", I4)')wann_ntype
write(100,*)
write(100,'("total number of WF : ", I4)')nwann
write(100,*)
do n=1,nwann
  iatom=iwann(1,n)
  lm=iwann(2,n)
  ispn=iwann(3,n)
  itype=iwann(4,n)
  write(100,'("  wf : ",I4)')n
  write(100,'("    type : ",I4)')itype
  write(100,'("    pure spinor orbital for projection : ")')
  write(100,'("      atom : ",I4)')iatom
  write(100,'("      l,m  : ",2I4)')lm2l(lm),lm-lm2l(lm)**2
  write(100,'("      ispn : ",I4)')ispn
  write(100,'("    energy interval : [",F8.4,",",F8.4,"]")')wann_eint(:,itype)
  write(100,'("    band interval : from ",I4," to ",I4)')wann_nint(:,itype)
  write(100,'("    potential : ",F8.4)')wann_v(itype)
  write(100,*)
enddo
close(100)

if (allocated(wann_c)) deallocate(wann_c)
allocate(wann_c(nwann,nstsv,nkptloc(iproc)))
wann_c=zzero

if (allocated(wann_unkmt)) deallocate(wann_unkmt)
allocate(wann_unkmt(lmmaxvr,nrfmax,natmtot,nspinor,nwann,nkptloc(iproc)))
wann_unkmt=zzero
if (allocated(wann_unkit)) deallocate(wann_unkit)
allocate(wann_unkit(ngkmax,nspinor,nwann,nkptloc(iproc)))
wann_unkit=zzero

if (allocated(wann_h)) deallocate(wann_h)
allocate(wann_h(nwann,nwann,nkpt))
wann_h=zzero
if (allocated(wann_e)) deallocate(wann_e)
allocate(wann_e(nwann,nkpt))
wann_e=0.d0
if (allocated(wann_p)) deallocate(wann_p)
allocate(wann_p(3,nwann,nwann,nkpt))
wann_p=zzero

if (allocated(wann_ene)) deallocate(wann_ene)
allocate(wann_ene(nwann))
wann_ene=0.d0
if (allocated(wann_occ)) deallocate(wann_occ)
allocate(wann_occ(nwann))
wann_occ=0.d0

if (allocated(wf_v_mtrx)) deallocate(wf_v_mtrx)
allocate(wf_v_mtrx(lmmaxlu,lmmaxlu,nspinor,nspinor,natmtot))

!if (allocated(wannmt)) deallocate(wannmt)
!allocate(wannmt(lmmaxvr,nrcmtmax,natmtot,nspinor,wann_nmax,nkptloc(iproc)))
!if (allocated(wannit)) deallocate(wannit)
!allocate(wannit(ngrtot,nspinor,wann_nmax,nkptloc(iproc)))

!if (allocated(wf_p)) deallocate(wf_p)
!allocate(wf_p(3,wf_dim,wf_dim,wann_nspin,nkpt))
!wf_p=dcmplx(0.d0,0.d0)

return
end
