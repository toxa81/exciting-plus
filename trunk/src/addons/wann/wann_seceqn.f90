subroutine wann_seceqn(ikloc,evecsv)
use modmain
use mod_wannier
! arguments
implicit none
integer, intent(in) :: ikloc
complex(8), intent(inout) :: evecsv(nstsv,nstsv)
! local variables
complex(8), allocatable :: z1(:,:),z2(:,:)
complex(8), allocatable :: work(:)
real(8), allocatable :: rwork(:)
integer i,j,i1,i2,info,lwork,itype,n,ik

ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)

lwork=2*nstsv
allocate(z1(nstsv,nstsv))
allocate(z2(nstsv,nstsv))
allocate(rwork(3*nstsv))
allocate(work(lwork))

z1=zzero
do i=1,nstsv
  z1(i,i)=evalsv(i,ik)
enddo

do n=1,nwantot
  itype=wan_info(4,n)
  do i1=1,nstsv
    do i2=1,nstsv
      z1(i1,i2)=z1(i1,i2)+dconjg(wann_c(n,i1,ikloc))*wann_c(n,i2,ikloc)*wann_v(itype)
    enddo
  enddo
enddo

if (ndmag.eq.1) then
! collinear: block diagonalise H
  call zheev('V','U',nstfv,z1(1,1),nstsv,evalsv(1,ik),work,lwork,rwork,info)
  if (info.ne.0) goto 20
  i=nstfv+1
  call zheev('V','U',nstfv,z1(i,i),nstsv,evalsv(i,ik),work,lwork,rwork,info)
  if (info.ne.0) goto 20
  do i=1,nstfv
    do j=1,nstfv
      z1(i,j+nstfv)=0.d0
      z1(i+nstfv,j)=0.d0
    end do
  end do
else
! non-collinear or spin-unpolarised: full diagonalisation
  call zheev('V','U',nstsv,z1,nstsv,evalsv(1,ik),work,lwork,rwork,info)
  if (info.ne.0) goto 20
endif
z2=zzero
do i1=1,nstsv
  do i2=1,nstsv
    do j=1,nstsv
      z2(i2,i1)=z2(i2,i1)+dconjg(z1(j,i1))*evecsv(i2,j)
    enddo
  enddo
enddo
evecsv=z2

deallocate(z1)
deallocate(z2)
deallocate(rwork)
deallocate(work)
return
20 continue
write(*,*)
write(*,'("Error(wann_seceqn): diagonalisation of the second-variational &
 &Hamiltonian failed")')
write(*,'(" for k-point ",I8)') ik
write(*,'(" ZHEEV returned INFO = ",I8)') info
write(*,*)
call pstop
return
end
