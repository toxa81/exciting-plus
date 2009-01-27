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

!if (allocated(wf_n)) deallocate(wf_n)
!allocate(wf_n(wf_dim,3))
!if (allocated(wf_lhbnd)) deallocate(wf_lhbnd)
!if (allocated(wf_lhen)) deallocate(wf_lhen)
!allocate(wf_lhbnd(2,wann_nspin,wf_dim))
!allocate(wf_lhen(2,wann_nspin,wf_dim))

!allocate(wf_deltav(wann_nspin,wf_dim))

do ispn=1,wann_nspin
  n=0
  do i=1,wann_natom
    iwgrp=wann_iatom(1+ispn,i)
    do lm=1,wann_iorbgrp(0,1,iwgrp)
      n=n+1
      iwann(n,ispn,1)=wann_iatom(1,i)
      iwann(n,ispn,2)=wann_iorbgrp(lm,1,iwgrp)
      iwann(n,ispn,3)=lm2l(iwann(n,ispn,2))
      iwann(n,ispn,4)=wann_iorbgrp(lm,2,iwgrp)
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
  enddo
enddo
close(100)
      
  
  
  
  
!  do lm=1,wann_iorb(0,i)
!    j=j+1
!    wf_n(j,1)=wann_iatom(i)
!    wf_n(j,2)=wann_iorb(lm,i)
!    wf_n(j,3)=lm2l(wann_iorb(lm,i))
!    do ispn=1,wann_nspin
!      wf_lhbnd(:,ispn,j)=wann_lhbnd(:,lm,ispn,i)
!      wf_lhen(:,ispn,j)=wann_lhen(:,lm,ispn,i)
!      wf_deltav(ispn,j)=wann_deltav(lm,ispn,i)   
!    enddo
!  enddo
!enddo


if (allocated(wann_c)) deallocate(wann_c)
allocate(wann_c(wann_nmax,nstfv,wann_nspin,nkptloc(iproc)))
wann_c=dcmplx(0.d0,0.d0)

if (allocated(wf_h)) deallocate(wf_h)
allocate(wf_h(wann_nmax,wann_nmax,wann_nspin,nkpt))
wf_h=dcmplx(0.d0,0.d0)

!if (allocated(wf_p)) deallocate(wf_p)
!allocate(wf_p(3,wf_dim,wf_dim,wann_nspin,nkpt))
!wf_p=dcmplx(0.d0,0.d0)

if (allocated(wf_e)) deallocate(wf_e)
allocate(wf_e(wann_nmax,wann_nspin,nkpt))
wf_e=0.d0

!if (allocated(wfpoco)) deallocate(wfpoco)
!allocate(wfpoco(nstsv,nstsv,nkptloc(iproc)))
!wfpoco=dcmplx(0.d0,0.d0)

!if (allocated(wfpoco1)) deallocate(wfpoco1)
!allocate(wfpoco1(nstsv,nstsv,nkptloc(iproc)))
!wfpoco1=dcmplx(0.d0,0.d0)

return
end
