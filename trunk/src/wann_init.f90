subroutine wann_init
use modmain
use modwann
implicit none

integer wann_natoms
integer, allocatable :: wann_iatom(:)
integer, allocatable :: wann_norb(:)
integer, allocatable :: wann_iorb(:,:)
integer, allocatable :: wann_lorb(:,:)
integer wann_maxnorb
integer i,j,lm

if (nspnfv.eq.2) then
  write(*,*)
  write(*,'("Error(wann_init): can''t do Wannier + spin spirals")')
  write(*,*)
  call pstop
endif

if (ndmag.eq.3) then
  write(*,*)
  write(*,'("Error(wann_init): can''t do Wannier + non-collinear magnetism")')
  write(*,*)
  call pstop
endif

wann_natoms=1
allocate(wann_iatom(wann_natoms))
wann_iatom(1)=1
allocate(wann_norb(wann_natoms))
wann_norb(1)=3
wann_maxnorb=maxval(wann_norb)
allocate(wann_iorb(wann_maxnorb,wann_natoms))
allocate(wann_lorb(wann_maxnorb,wann_natoms))
wann_iorb(1:3,1)=(/5,6,8/)
wann_lorb(1:3,1)=(/2,2,2/)

wf_dim=0
do i=1,wann_natoms
  wf_dim=wf_dim+wann_norb(i)
enddo
if (allocated(wf_n)) deallocate(wf_n)
allocate(wf_n(wf_dim,6))
if (allocated(wf_lhbnd)) deallocate(wf_lhbnd)
allocate(wf_lhbnd(wf_dim,2))
j=0
do i=1,wann_natoms
  do lm=1,wann_norb(i)
    j=j+1
    wf_n(j,1)=wann_iatom(i)
    wf_n(j,2)=wann_iorb(lm,i)
    wf_n(j,3)=wann_lorb(lm,i)
  enddo
enddo

open(100,file='WANNIER.OUT',form='formatted',status='replace')
write(100,*)'wf_dim=',wf_dim
do i=1,wf_dim
  write(100,*)'wf=',i,' atom=',wf_n(i,1),' lm=',wf_n(i,2),' l=',wf_n(i,3)
enddo
close(100)

wf_lhbnd(1:wf_dim,1)=8
wf_lhbnd(1:wf_dim,2)=12

if (allocated(a_ort)) deallocate(a_ort)
allocate(a_ort(wf_dim,nstfv,nspinor,nkpt))
a_ort=dcmplx(0.d0,0.d0)

if (allocated(wf_h)) deallocate(wf_h)
allocate(wf_h(wf_dim,wf_dim,nspinor,nkpt))
wf_h=dcmplx(0.d0,0.d0)

if (allocated(wf_e)) deallocate(wf_e)
allocate(wf_e(wf_dim,nspinor,nkpt))
wf_e=dcmplx(0.d0,0.d0)

return
end
