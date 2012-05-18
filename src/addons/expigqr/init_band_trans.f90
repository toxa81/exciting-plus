subroutine init_band_trans
use modmain
use mod_wannier
use mod_expigqr
use mod_nrkp
implicit none
integer ikloc,i,ik,jk,ist1,ist2
logical, external :: wann_diel
!
if (allocated(nmegqblh)) deallocate(nmegqblh)
allocate(nmegqblh(nkptnrloc))
nmegqblh=0
if (allocated(bmegqblh)) deallocate(bmegqblh)
allocate(bmegqblh(2,nstsv*nstsv,nkptnrloc))
bmegqblh=0
if (wannier_megq) then
  if (allocated(nmegqblhwan)) deallocate(nmegqblhwan)
  allocate(nmegqblhwan(nkptnrloc))
  nmegqblhwan=0
  if (allocated(imegqblhwan)) deallocate(imegqblhwan)
  allocate(imegqblhwan(nstsv*nstsv,nkptnrloc))
  imegqblhwan=0
endif
call getmeidx
if (wannier_megq) then
  nmegqblhwanmax=maxval(nmegqblhwan)
endif
! get minimal transition energy
lr_min_e12=1.d10
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  jk=idxkq(1,ik)
  do i=1,nmegqblh(ikloc)
    ist1=bmegqblh(1,i,ikloc)
    ist2=bmegqblh(2,i,ikloc)
    lr_min_e12=min(lr_min_e12,abs(evalsvnr(ist1,ik)-evalsvnr(ist2,jk)))
  enddo
enddo
call mpi_grid_reduce(lr_min_e12,dims=(/dim_k/),op=op_min)
return
end
