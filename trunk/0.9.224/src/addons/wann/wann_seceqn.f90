subroutine wann_seceqn(ik,evecsv)
use modmain
! arguments
implicit none
integer, intent(in) :: ik
complex(8), intent(inout) :: evecsv(nstsv,nstsv)
! local variables
complex(8), allocatable :: z1(:,:),z2(:,:)
real(8), allocatable :: work(:)
complex(8), allocatable :: rwork(:)
integer i,j,n,ispn,itype,i1,i2,info,lwork
integer, external :: ikglob

lwork=2*nstsv
allocate(z1(nstsv,nstsv))
allocate(z2(nstsv,nstsv))
allocate(rwork(3*nstsv))
allocate(work(lwork))

z1=zzero
do i=1,nstsv
  z1(i,i)=evalsv(i,ikglob(ik))
enddo

do ispn=1,nspinor
  do n=1,nwann(ispn)
    itype=iwann(n,ispn,4)
    do i1=1,nstfv
      do i2=1,nstfv
        z1(i1+(ispn-1)*nstfv,i2+(ispn-1)*nstfv)=z1(i1+(ispn-1)*nstfv,i2+(ispn-1)*nstfv)+&
          dconjg(wann_c(n,i1,ispn,ik))*wann_c(n,i2,ispn,ik)*wann_v(itype)
      enddo
    enddo
  enddo
enddo

! collinear: block diagonalise H
call zheev('V','U',nstfv,z1,nstsv,evalsv(1,ikglob(ik)),work,lwork,rwork,info)
if (info.ne.0) goto 20
i=nstfv+1
call zheev('V','U',nstfv,z1(i,i),nstsv,evalsv(i,ikglob(ik)),work,lwork,rwork,info)
if (info.ne.0) goto 20
do i=1,nstfv
  do j=1,nstfv
    z1(i,j+nstfv)=0.d0
    z1(i+nstfv,j)=0.d0
  end do
end do

z2=zzero
do i1=1,nstsv
  do i2=1,nstsv
    do j=1,nstsv
      z2(i2,i1)=z2(i2,i1)+dconjg(z1(j,i1))*evecsv(i2,j)
    enddo
  enddo
enddo
evecsv=z2

deallocate(z1,z2,rwork,work)
return


20 continue
write(*,*)
write(*,'("Error(wann_seceqn): diagonalisation of the second-variational &
 &Hamiltonian failed")')
write(*,'(" for k-point ",I8)') ikglob(ik)
write(*,'(" ZHEEV returned INFO = ",I8)') info
write(*,*)
stop


return
end