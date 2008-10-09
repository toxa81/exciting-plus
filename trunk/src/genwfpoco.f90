subroutine genwfpoco(ik)
use modmain
use modwann
implicit none
! arguments
integer, intent(in) :: ik

integer i,j,n,ispn,istfv
complex(8) zt1

complex(8), allocatable :: wtmp(:,:,:)

! calculate <u_{n,k,\sigma}|\phi_{i}>, where |u_{n,k,\sigma}> are Fourier-transforms of 
!  Wannier functions and |\phi_{i}> are first-variational states
allocate(wtmp(wf_dim,wann_nspins,nstfv))
wtmp=dcmplx(0.d0,0.d0)
do n=1,wf_dim
  do ispn=1,wann_nspins
    do istfv=1,nstfv
      do j=1,nstfv
        wtmp(n,ispn,istfv)=wtmp(n,ispn,istfv)+wfc(n,j,ispn,ik) * &
	  evecsvloc(istfv+(ispn-1)*nstfv,j+(ispn-1)*nstfv,ik)
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

return
end
