subroutine writenn
use modmain
implicit none

integer              :: i, iat, is, js, jat, i1, i2, i3, tmp(6)

integer              :: nneigh
integer ,allocatable :: clust(:,:)

!--- radius of cluster of nearest neighbours (in a.u.)
real*8               :: rclust = 10.d0
real*8               :: a(3),v1(3),d1,v2(3)

integer              :: llim(3)

open(50,file='NN.OUT',status='replace',form='formatted')

do i = 1, 3
  a(i) = sqrt(avec(1,i)**2+avec(2,i)**2+avec(3,i)**2)
  llim(i) = rclust/a(i) + 1
  write(50,'("avec(",I1") = ",3G18.10," length=",G18.10," llim=",I2)')i,avec(:,i),a(i),llim(i)
enddo

do is = 1, nspecies
  do iat = 1, natoms(is)
    write(50,'("is=",I2," iat=",I2," lat.coord.=",3G18.10," cart.coord.=",3G18.10)') &
      is,iat,atposl(:,iat,is),atposc(:,iat,is)
  enddo
enddo

allocate(clust(natmtot*(2*llim(1)+1)*(2*llim(2)+1)*(2*llim(3)+1),6))

do is = 1, nspecies
  do iat = 1, natoms(is)
    clust = 0
    nneigh = 0
    
    do js = 1, nspecies 
      do jat = 1, natoms(js)
        do i1 = -llim(1), llim(1)
	do i2 = -llim(2), llim(2)
	do i3 = -llim(3), llim(3)
	  v1(:) = atposc(:,jat,js) + i1*avec(:,1)+i2*avec(:,2)+i3*avec(:,3) - atposc(:,iat,is)
	  d1 = sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
	  if (d1.le.rclust) then
	    nneigh = nneigh + 1
	    clust(nneigh,1) = js
	    clust(nneigh,2) = jat
	    clust(nneigh,3) = i1
	    clust(nneigh,4) = i2
	    clust(nneigh,5) = i3
	    clust(nneigh,6) = d1*100000
	  endif
	enddo
	enddo
	enddo
      enddo
    enddo
    
    do i1 = 1, nneigh - 1
    do i2 = i1+1, nneigh
      if (clust(i1,6).gt.clust(i2,6)) then
        tmp(:) = clust(i1,:)
	clust(i1,:) = clust(i2,:)
	clust(i2,:) = tmp(:)
      endif
    enddo
    enddo
    
    write(50,'("Cluster arount iat=",I2,"(is=",A,")")')iat,trim(spsymb(is))
    do i = 1, nneigh
      js = clust(i,1)
      jat = clust(i,2)
      v1 = atposl(:,jat,js) + clust(i,3:5)
      v2 = atposc(:,jat,js) + clust(i,3)*avec(:,1)+clust(i,4)*avec(:,2)+clust(i,5)*avec(:,3) - atposc(:,iat,is)

      write(50,'(I4," jat=",I2,"(",A2,") T=",3I3," dist=",F12.6," pos(lat)=",3F12.6," pos(Cart)=",3F12.6)') &
        i,clust(i,2),trim(spsymb(clust(i,1))),clust(i,3:5),clust(i,6)/100000.d0,v1,v2
    enddo
  enddo
enddo    

close(50)

deallocate(clust)


end
