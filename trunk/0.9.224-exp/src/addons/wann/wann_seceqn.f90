subroutine wann_seceqn(ik,evecsv)
use modmain
! arguments
implicit none
integer, intent(in) :: ik
complex(8), intent(inout) :: evecsv(nstsv,nstsv)
! local variables
complex(8), allocatable :: z1(:,:),z2(:,:)
complex(8), allocatable :: work(:)
real(8), allocatable :: rwork(:)
integer i,j,n,ispn,itype,i1,i2,info,lwork,n1,n2,ias,l,m1,m2,lm1,lm2
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

do n=1,nwann
  itype=iwann(4,n)
  do i1=1,nstsv
    do i2=1,nstsv
      z1(i1,i2)=z1(i1,i2)+dconjg(wann_c(n,i1,ik))*wann_c(n,i2,ik)*wann_v(itype)
    enddo
  enddo
enddo

!do i=1,wann_natom
!  ias=wann_iprj(1,i)
!  do l=0,lmaxlu
!    do m1=-l,l
!      do m2=-l,l
!        lm1=idxlm(l,m1)
!        lm2=idxlm(l,m2)
!        n1=iasiwann(ias,lm1,ispn)
!        n2=iasiwann(ias,lm2,ispn)
!        if (n1.ne.-1.and.n2.ne.-1) then
!          do i1=1,nstfv
!            do i2=1,nstfv
!              z1(i1+(ispn-1)*nstfv,i2+(ispn-1)*nstfv)=z1(i1+(ispn-1)*nstfv,i2+(ispn-1)*nstfv)+&
!                dconjg(wann_c(n1,i1,ispn,ik))*wann_c(n2,i2,ispn,ik)*wf_v_mtrx(lm1,lm2,ispn,ispn,ias)
!            enddo
!          enddo
!        endif
!      enddo !m2
!    enddo !m1
!  enddo !l
!enddo !i

if (ndmag.eq.1) then
  ! collinear: block diagonalise H
  call zheev('V','U',nstfv,z1(1,1),nstsv,evalsv(1,ikglob(ik)),work,lwork,rwork,info)
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
else
! non-collinear or spin-unpolarised: full diagonalisation
  call zheev('V','U',nstsv,z1,nstsv,evalsv(1,ikglob(ik)),work,lwork,rwork,info)
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
write(*,'(" for k-point ",I8)') ikglob(ik)
write(*,'(" ZHEEV returned INFO = ",I8)') info
write(*,*)
call pstop
return
end