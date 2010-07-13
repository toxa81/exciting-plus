subroutine genwann_h(ikloc)
use modmain
use mod_mpi_grid
implicit none
integer, intent(in) :: ikloc
complex(8), allocatable :: zt2(:,:)
integer m1,m2,j
integer ik
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)

! compute H(k) 
allocate(zt2(nwann,nwann))
zt2=zzero
if (nwann_h.eq.0) then
  do m1=1,nwann
    do m2=1,nwann
      do j=1,nstsv
        zt2(m1,m2)=zt2(m1,m2)+dconjg(wann_c(m1,j,ikloc))*wann_c(m2,j,ikloc) * &
          evalsv(j,ik)
      enddo
    enddo
  enddo
else
  do m1=1,nwann_h
    do m2=1,nwann_h
      do j=1,nstsv
        zt2(m1,m2)=zt2(m1,m2)+dconjg(wann_c(iwann_h(m1),j,ikloc))*&
          wann_c(iwann_h(m2),j,ikloc)*evalsv(j,ik)
      enddo
    enddo
  enddo
  do m1=nwann_h+1,nwann
    zt2(m1,m1)=-100.d0*zone
  enddo
endif
wann_h(:,:,ik)=zt2(:,:)
call diagzhe(nwann,zt2,wann_e(1,ik))
deallocate(zt2)
return
end
