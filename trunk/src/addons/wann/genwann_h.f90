subroutine genwann_h(ikloc)
use modmain
implicit none
integer, intent(in) :: ikloc
complex(8), allocatable :: zt2(:,:)
integer m1,m2,j
integer, external :: ikglob

! compute H(k) 
allocate(zt2(nwann,nwann))
zt2=zzero
do m1=1,nwann
  do m2=1,nwann
    do j=1,nstsv
      zt2(m1,m2)=zt2(m1,m2)+dconjg(wann_c(m1,j,ikloc))*wann_c(m2,j,ikloc) * &
        evalsv(j,ikglob(ikloc))
    enddo
  enddo
enddo
wann_h(:,:,ikglob(ikloc))=zt2(:,:)
call diagzhe(nwann,zt2,wann_e(1,ikglob(ikloc)))
deallocate(zt2)
return
end
