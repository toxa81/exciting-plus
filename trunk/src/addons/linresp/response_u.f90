subroutine response_u
use modmain
implicit none

complex(8), allocatable :: uscrn(:,:)
complex(8), allocatable :: ubare(:,:)
complex(8), allocatable :: z1(:,:),z2(:,:)
integer iq,i,j
character*100 fname,qnm

allocate(uscrn(nwann,nwann))
allocate(ubare(nwann,nwann))
allocate(z1(nwann,nwann))
allocate(z2(nwann,nwann))

uscrn=zzero
ubare=zzero

do iq=1,nvq0
  call qname(ivq0m_list(:,iq),qnm)
  qnm="./"//trim(qnm)//"/"//trim(qnm)
  fname=trim(qnm)//"_U"
  open(170,file=trim(fname),status='old',form='unformatted')
  read(170)z1,z2
  close(170)
  uscrn=uscrn+z1
  ubare=ubare+z2
enddo

uscrn=uscrn/nvq0
ubare=ubare/nvq0

open(150,file='CRPA_U.OUT',status='replace',form='formatted')
write(150,'("Screened U matrix")')
write(150,'("real part")')
do i=1,nwann
  write(150,'(100F12.6)')(dreal(uscrn(i,j)),j=1,nwann)
enddo
write(150,'("imag part")')
do i=1,nwann
  write(150,'(100F12.6)')(dimag(uscrn(i,j)),j=1,nwann)
enddo
write(150,*)
write(150,'("Bare U matrix")')
write(150,'("real part")')
do i=1,nwann
  write(150,'(100F12.6)')(dreal(ubare(i,j)),j=1,nwann)
enddo
write(150,'("imag part")')
do i=1,nwann
  write(150,'(100F12.6)')(dimag(ubare(i,j)),j=1,nwann)
enddo
close(150)
deallocate(uscrn,ubare,z1,z2)

return
end


