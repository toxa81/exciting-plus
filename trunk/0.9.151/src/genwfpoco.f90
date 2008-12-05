subroutine genwfpoco(ik,evecsv)
use modmain
use modwann
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: evecsv(nstsv,nstsv)

integer i,j,n,ispn,istfv
complex(8) zt1

complex(8), allocatable :: wtmp(:,:,:)

wfpoco1(:,:,ik)=wfpoco(:,:,ik)
! calculate <u_{n,k,\sigma}|\phi_{i}>, where |u_{n,k,\sigma}> are Fourier-transforms of 
!  Wannier functions and |\phi_{i}> are first-variational states
allocate(wtmp(wf_dim,wann_nspins,nstfv))
wtmp=dcmplx(0.d0,0.d0)
do n=1,wf_dim
  do ispn=1,wann_nspins
    do istfv=1,nstfv
      do j=1,nstfv
        wtmp(n,ispn,istfv)=wtmp(n,ispn,istfv)+wfc(n,j,ispn,ik) * &
	  evecsv(istfv+(ispn-1)*nstfv,j+(ispn-1)*nstfv)
      enddo
    enddo
  enddo
enddo
do i=1,nstfv
  do j=1,nstfv
    do ispn=1,wann_nspins
      zt1=dcmplx(0.d0,0.d0)
      do n=1,wf_dim
        zt1=zt1+wtmp(n,ispn,i)*dconjg(wtmp(n,ispn,j))*wf_deltav(ispn,n)
      enddo
      wfpoco(i+(ispn-1)*nstfv,j+(ispn-1)*nstfv,ik)=zt1
    enddo
  enddo
enddo
deallocate(wtmp)

!wfpoco(:,:,ik)=0.1d0*wfpoco(:,:,ik)+0.9d0*wfpoco1(:,:,ik)

return
end

subroutine writewfpoco
use modmain
use modwann
implicit none

integer ik,i,j1,j2
character*4 c4
integer, external :: ikglob

if (iproc.eq.0) then
  call write_integer(nkpt,1,'/dimensions','nkpt')
  call write_integer(nstsv,1,'/dimensions','nstsv')
endif
  
do i=0,nproc-1
  if (i.eq.iproc) then
    do ik=1,nkptloc(iproc)
      write(c4,'(I4.4)')ikglob(ik)
      call write_complex16(wfpoco(1,1,ik),nstsv*nstsv,'/kpoints/'//c4,'wfpoco')
      call write_real8(vkl(1,ikglob(ik)),3,'/kpoints/'//c4,'vkl')
    enddo
  endif
  call barrier
enddo
return
end

subroutine readwfpoco
use modmain
use modwann
implicit none

integer ik,i,j1,j2,nkpt_,nstsv_
integer, external :: ikglob
character*4 c4

do i=0,nproc-1
  if (i.eq.iproc) then
    do ik=1,nkptloc(iproc)
      write(c4,'(I4.4)')ikglob(ik)
      call read_complex16(wfpoco(1,1,ik),nstsv*nstsv,'/kpoints/'//c4,'wfpoco')
    enddo
  endif
  call barrier
enddo
return
end



