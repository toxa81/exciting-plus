subroutine writenn
use modmain
implicit none

integer i,ia,is,ias,ja,js,jas,i1,i2,i3,tmp(6)

integer nneigh
integer ,allocatable :: clust(:,:)

! radius of cluster of nearest neighbours (in a.u.)
real(8) rclust

real(8) a(3),v1(3),d1,v2(3),v3(3),dst
integer llim(3)

rclust=20.d0

open(50,file='NN.OUT',status='replace',form='formatted')

write(50,*)
write(50,'("Radius of cluster of nearest neighbours : ",G18.10)')rclust

write(50,*)
write(50,'("Lattice vectors and limits : ")')
! find lattice limits
do i=1,3
  a(i)=sqrt(avec(1,i)**2+avec(2,i)**2+avec(3,i)**2)
  llim(i)=rclust/a(i)+1
  write(50,'("avec(",I1")=",3G18.10," length=",G18.10," llim=",I2)')i,avec(:,i),a(i),llim(i)
enddo

write(50,*)
write(50,'("Atom list : ")')
write(50,'(" iatom (is ia ias)",15(" "),"pos(lat)",24(" "),"pos(Cart)")')
write(50,'(85("-"))')
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,'(3X,A,T6," (",I2,1X,I2,1X,I3,")",T22,3F10.5,2X,3F10.5)') &
      trim(spsymb(is)),is,ia,ias,atposl(:,ia,is),atposc(:,ia,is)
  enddo
enddo

allocate(clust(natmtot*(2*llim(1)+1)*(2*llim(2)+1)*(2*llim(3)+1),6))

do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    
    clust=0
    nneigh=0
    do js=1,nspecies 
      do ja=1,natoms(js)
        do i1=-llim(1),llim(1)
	do i2=-llim(2),llim(2)
	do i3=-llim(3),llim(3)
	  v1(:)=atposc(:,ja,js)+i1*avec(:,1)+i2*avec(:,2)+i3*avec(:,3)-atposc(:,ia,is)
	  d1=sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
	  if (d1.le.rclust) then
	    nneigh=nneigh+1
	    clust(nneigh,1)=js
	    clust(nneigh,2)=ja
	    clust(nneigh,3)=i1
	    clust(nneigh,4)=i2
	    clust(nneigh,5)=i3
	    clust(nneigh,6)=d1*100000
	  endif
	enddo
	enddo
	enddo
      enddo
    enddo
    
    do i1=1,nneigh-1
    do i2=i1+1,nneigh
      if (clust(i1,6).gt.clust(i2,6)) then
        tmp(:)=clust(i1,:)
	clust(i1,:)=clust(i2,:)
	clust(i2,:)=tmp(:)
      endif
    enddo
    enddo
    
    write(50,*)
    write(50,'("Cluster around ",A," (is ia ias : ",I2,1X,I2,1X,I3")")')trim(spsymb(is)),is,ia,ias
    write(50,'(" jatom (js ja jas)   D(a.u.)     D(A)",8(" "),"T",18(" "),"R(Cart)")')
    write(50,'(82("-"))')
    do i=1,nneigh
      js=clust(i,1)
      ja=clust(i,2)
      jas=idxas(ja,js)
      v1=atposl(:,ja,js)+clust(i,3:5)
      v2=atposc(:,ja,js)+clust(i,3)*avec(:,1)+clust(i,4)*avec(:,2)+clust(i,5)*avec(:,3)
      v3=v2(:)-atposc(:,ia,is)
      dst=clust(i,6)/100000.d0
      write(50,'(3X,A,T6," (",I2,1X,I2,1X,I3") ",2F10.5,2X,3I3,2X,3F10.5)') &
        trim(spsymb(clust(i,1))),js,ja,jas,dst,dst*au2ang,clust(i,3:5),v3
    enddo
  enddo
enddo    

close(50)

deallocate(clust)


end
