subroutine wann_init
use modmain
use modwann
implicit none

integer i,j,lm,ispn

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

wf_dim=0
do i=1,wann_natoms
  wf_dim=wf_dim+wann_iorb(0,i)
enddo
if (allocated(wf_n)) deallocate(wf_n)
allocate(wf_n(wf_dim,3))
if (allocated(wf_lhbnd)) deallocate(wf_lhbnd)
if (allocated(wf_lhen)) deallocate(wf_lhen)
allocate(wf_lhbnd(2,wann_nspins,wf_dim))
allocate(wf_lhen(2,wann_nspins,wf_dim))

allocate(wf_deltav(wann_nspins,wf_dim))

j=0
do i=1,wann_natoms
  do lm=1,wann_iorb(0,i)
    j=j+1
    wf_n(j,1)=wann_iatom(i)
    wf_n(j,2)=wann_iorb(lm,i)
    wf_n(j,3)=lm2l(wann_iorb(lm,i))
    do ispn=1,wann_nspins
      wf_lhbnd(:,ispn,j)=wann_lhbnd(:,lm,ispn,i)
      wf_lhen(:,ispn,j)=wann_lhen(:,lm,ispn,i)
      wf_deltav(ispn,j)=wann_deltav(lm,ispn,i)   
    enddo
  enddo
enddo

open(100,file='WANNIER.OUT',form='formatted',status='replace')
write(100,*)'wf_dim=',wf_dim
do i=1,wf_dim
  write(100,*)'wf=',i,' atom=',wf_n(i,1),' lm=',wf_n(i,2),' l=',wf_n(i,3),'lh_bnd=',wf_lhbnd(:,1,i),'deltav=',wf_deltav(1,i)
enddo
close(100)

if (allocated(wfc)) deallocate(wfc)
allocate(wfc(wf_dim,nstfv,wann_nspins,nkptloc(iproc)))
wfc=dcmplx(0.d0,0.d0)

if (allocated(wf_h)) deallocate(wf_h)
allocate(wf_h(wf_dim,wf_dim,wann_nspins,nkpt))
wf_h=dcmplx(0.d0,0.d0)

if (allocated(wf_e)) deallocate(wf_e)
allocate(wf_e(wf_dim,wann_nspins,nkpt))
wf_e=0.d0

if (allocated(wfpoco)) deallocate(wfpoco)
allocate(wfpoco(nstsv,nstsv,nkptloc(iproc)))
wfpoco=dcmplx(0.d0,0.d0)

return
end
