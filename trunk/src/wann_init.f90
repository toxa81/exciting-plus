subroutine wann_init
use modmain
use modwann
implicit none

if (nspnfv.eq.2) then
  write(*,*)
  write(*,'("Error(wann_init): can''t do Wannier + spin spirals")')
  write(*,*)
  call pstop
endif

wannier=.true.

wf_dim=8
if (allocated(wf_n)) deallocate(wf_n)
allocate(wf_n(wf_dim,6))
if (allocated(wf_lhbnd)) deallocate(wf_lhbnd)
allocate(wf_lhbnd(wf_dim,2))

wf_n(1:5,1)=1
wf_n(6:8,1)=2
wf_n(1,2)=5
wf_n(2,2)=6
wf_n(3,2)=7
wf_n(4,2)=8
wf_n(5,2)=9
wf_n(6,2)=2
wf_n(7,2)=3
wf_n(8,2)=4

wf_n(1:5,3)=2
wf_n(6:8,3)=1

wf_lhbnd(1:wf_dim,1)=5
wf_lhbnd(1:wf_dim,2)=12

if (allocated(a_ort)) deallocate(a_ort)
allocate(a_ort(wf_dim,nstsv,nspinor,nkpt))
a_ort=dcmplx(0.d0,0.d0)

if (allocated(wf_h)) deallocate(wf_h)
allocate(wf_h(wf_dim,wf_dim,nspinor,nkpt))
wf_h=dcmplx(0.d0,0.d0)

if (allocated(wf_e)) deallocate(wf_e)
allocate(wf_e(wf_dim,nspinor,nkpt))
wf_e=dcmplx(0.d0,0.d0)

return
end
