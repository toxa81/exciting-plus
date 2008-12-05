subroutine response_jdos
use modmain
implicit none
real(8), allocatable :: occsvnr(:,:)
real(8), allocatable :: evalsvnr(:,:)
! number of energy-mesh points
integer nepts
! energy mesh
complex(8), allocatable :: w(:)
integer ik,i,n1,n2
complex(8), allocatable :: f(:)

allocate(occsvnr(nstsv,nkptnr))
allocate(evalsvnr(nstsv,nkptnr))


do ik=1,nkptnr
  call getoccsv(vklnr(1,ik),occsvnr(1,ik))
  call getevalsv(vklnr(1,ik),evalsvnr(1,ik))
enddo

! setup energy mesh
nepts=1+maxomega/domega
allocate(w(nepts))
do i=1,nepts
  w(i)=dcmplx(domega*(i-1),eta)/ha2ev
enddo

allocate(f(nepts))
f=dcmplx(0.d0,0.d0)
do i=1,nepts
  do n1=1,nstfv
    do n2=1,nstfv
      do ik=1,nkptnr
        f(i)=f(i)+(occsvnr(n1,ik)-occsvnr(n2,ik))/(evalsvnr(n1,ik)-evalsvnr(n2,ik)+w(i))
      enddo
    enddo
  enddo
enddo
f=f/nkptnr/ha2ev

open(50,file='JDOS.OUT',form='formatted',status='replace')
do i=1,nepts
  write(50,'(2F12.6)')dreal(w(i))*ha2ev,-dimag(f(i))
enddo
close(50)
  
return
end
