subroutine writebz
use modmain
implicit none
real(8) pt1(3,26),p0(3),pt2(3,26)
integer i1,i2,i3,i,m,npt,nc,nf,n1
integer connections(2,100)
integer faces(0:10,100)
logical vbz,l1

call readinput
call init0

! connect center point with neigbours
m=1
do i1=-1,1
  do i2=-1,1
    do i3=-1,1
      if (.not.(i1.eq.0.and.i2.eq.0.and.i3.eq.0)) then
        pt1(:,m)=i1*bvec(1,:)+i2*bvec(2,:)+i3*bvec(3,:)
        m=m+1
    endif
    enddo
  enddo
enddo

! find intersection of the planes, perpendicular to pt1 and
!   going through midpoints; this will define the points of BZ polyhedra
npt=0
do i1=1,26
  do i2=1,26
    do i3=1,26 
      call find_intersection(pt1(:,i1)/2,pt1(:,i2)/2,pt1(:,i3)/2,p0)
      call is_vector_in_bz(p0,pt1,vbz)
      if (vbz) then
! check if the point is already added to the list      
        l1=.true.
        do i=1,npt
          if (sum(abs(p0(:)-pt2(:,i))).lt.1d-6) l1=.false.
        enddo
        if (l1) then
          npt=npt+1
          pt2(:,npt)=p0(:)
        endif
      endif
    enddo
  enddo
enddo
call connect_pts(pt2,npt,connections,nc,faces,nf)
! count points starting from 0
faces(1:10,:)=faces(1:10,:)-1
connections(:,:)=connections(:,:)-1

n1=0
do i=1,nf
  n1=n1+faces(0,i)
enddo

open(150,file='bz.dx',status='replace',form='formatted')

write(150,'("object 1 class array type string rank 1 shape 5 items",\
  1x,I3,1x,"data follows")')npt
do i=1,npt
  write(150,'(""",I2,""")')i-1
enddo
close(150)
!write(64,14)

!write(64,16)np
!do i = 1, np
!  write(64,18)0.d0,0.d0,1.d0
!enddo
!write(64,14)
!
!write(64,20)np
!do i = 1, np
!  write(64,18)bz_pt(i,:)
!enddo
!write(64,14)
!
!write(64,22)nc
!do i = 1, nc
!  write(64,24)connections(i,:)
!enddo
!write(64,26)
!write(64,28)
!
!write(64,30)n1
!do i = 1, nf
!  write(64,32)faces(i,1:faces(i,0))
!enddo
!write(64,28)
!
!write(64,34)nf
!i = 0
!j = 1
!do while ( i.lt.n1 )
!  write(64,36)i
!i = i + faces(j,0)
!j = j + 1
!enddo
!write(64,38)
!
!write(64,40)nf
!do i = 1, nf
!  write(64,36)i-1
!enddo
!write(64,42)
!
!write(64,60)
!
!close(64)


10  format('object 1 class array type string rank 1 shape 5 items',1x,i3,1x,'data follows')
12  format('"',i2,'"')
14  format('attribute "dep" string "positions"')
16  format('object 2 class array type float rank 1 shape 3 items',1x,i3,1x,'data follows')
18  format(3f12.6)
20  format('object 3 class array type float rank 1 shape 3 items',1x,i3,1x,'data follows')
22  format('object 4 class array type int rank 1 shape 2 items',1x,i3,1x,'data follows')
24  format(2i3)
26  format('attribute "element type" string "lines"')
28  format('attribute "ref" string "positions"')
30  format('object 5 class array type int rank 0 items',1x,i3,1x,'data follows')
32  format(10i3)
34  format('object 6 class array type int rank 0 items',1x,i3,1x,'data follows')
36  format(i3)
38  format('attribute "ref" string "edges"')
40  format('object 7 class array type int rank 0 items',1x,i3,1x,'data follows')
42  format('attribute "ref" string "loops"')
60  format('object "bz" class field'          &
       /'component "data" value 1'        &
   /'component "colors" value 2'      &
   /'component "positions" value 3'   &
   /'component "connections" value 4' &
   /'component "edges" value 5'       & 
   /'component "loops" value 6'       &
   /'component "faces" value 7'       &
   /'end')

end

return
end


subroutine find_intersection(q1,q2,q3,p0)
implicit none
! arguments
real(8), intent(in) :: q1(3)
real(8), intent(in) :: q2(3)
real(8), intent(in) :: q3(3)
real(8), intent(out) :: p0(3)
! local variables
real(8) :: a(3,3),b(3),a1(3,3)
real(8) :: det0,det1,det2,det3
real(8), external :: r3mdet

! equation for a plane: \vec{n}*(\vec{x}-\vec{q})=0
!
! we are looking for the point p0 which satisfies three equations
!
! for the particular task of BZ construction normal of the plane is 
!   parallel to the direction from center point to nearest neighbour 
!   lattice point, so we may put \vec{n}=\vec{q} 

! a-matrix for the system of linear equations
a(1,:)=q1(:)
a(2,:)=q2(:)
a(3,:)=q3(:)

det0=r3mdet(a)
! linear dependent normals -> intersection of 3 planes is not a point
if (abs(det0).lt.1d-6) then
  p0(:)=1000.d0
  return
endif

! b-vector for the system of linear equations
b(1)=(q1(1)**2+q1(2)**2+q1(3)**2)
b(2)=(q2(1)**2+q2(2)**2+q2(3)**2)
b(3)=(q3(1)**2+q3(2)**2+q3(3)**2)

! solve system of linear equations using Cramer's rule
a1(:,:)=a(:,:)
a1(:,1)=b(:)
det1=r3mdet(a1)

a1(:,:)=a(:,:)
a1(:,2)=b(:)
det2=r3mdet(a1)

a1(:,:)=a(:,:)
a1(:,3)=b(:)
det3=r3mdet(a1)

p0(:)=(/det1/det0,det2/det0,det3/det0/)

return
end


subroutine is_vector_in_bz(p0,pt1,f)
implicit none
real(8), intent(in) :: p0(3)
real(8), intent(in) :: pt1(3,26)
logical, intent(out) :: f

real(8) q1,q2
integer i

! distance from the center of BZ to the point p0
q1=sqrt(p0(1)**2+p0(2)**2+p0(3)**2)-1.d-4
f=.true.

do i=1,26
! distance from the nearest neigbour lattice point to the point p0
  q2=sqrt((p0(1)-pt1(1,i))**2+(p0(2)-pt1(2,i))**2+(p0(3)-pt1(3,i))**2)
! check the main property of BZ: for any point inside BZ distance to center 
!   is smaller than distance to neighbouring lattice points
  if (q2.lt.q1) f=.false.
enddo

return
end


subroutine connect_pts(points,np,connections,nc,faces,nf)
implicit   none
! arguments
integer, intent(in) :: np
real(8), intent(in) :: points(3,np)
integer, intent(out) :: nc
integer, intent(out) :: connections(2,100)
integer, intent(out) :: nf
integer, intent(out) :: faces(0:10,100)

integer i1,i2,i3,m,n,np1,j,j1
real(8) pts_in_plane(3,10),norm(3),q1(3),q2(3),q3(3),q4(3)
integer pts_idx(10),face(0:10),conn(2,10)
real(8) fpt(np),a1,a2,a3
logical l1

nc=0
nf=0

do i1=1,np-2
  do i2=i1+1,np-1
    do i3=i2+1,np
! make a plane from any three poins and calculate plane's normal
      q1(:)=points(:,i2)-points(:,i1)
      q2(:)=points(:,i3)-points(:,i1)
! norm is a vector product
      norm(1)=q1(2)*q2(3)-q2(2)*q1(3)
      norm(2)=-(q1(1)*q2(3)-q2(1)*q1(3)) 
      norm(3)=q1(1)*q2(2)-q2(1)*q1(2)
! calculate scalar product between normal and vector q3
!  vector q3 is the difference between any point (index m) and point 
!  belonging to a plane (index i1)
      do m=1,np
        q3(:)=points(:,m)-points(:,i1)
        fpt(m)=(norm(1)*q3(1)+norm(2)*q3(2)+norm(3)*q3(3))
      enddo
! if all scalar products have the same sign, then all poins are on the same side
!  of the plane -> points i1,i2,i3 belong to the face of polyhedra 
      if (all(fpt(:).gt.-1d-6).or.all(fpt(:).lt.1d-6)) then
! norm should point out of BZ (to the side of the plane without points)
        if (all(fpt(:).gt.-1d-6)) norm(:)=-norm(:)
! find all points lying in the plane; this are the points, satisfying the plane
!  equation n*(x-p)=0
        np1=0
        do m=1,np
          if (abs(fpt(m)).lt.1d-6) then
            np1=np1+1
            pts_in_plane(:,np1)=points(:,m)
            pts_idx(np1)=m
          endif
        enddo !m
! find the average point (center of the face)
        q3(:)=0.d0
        do m=1,np1
          q3(:)=q3(:)+pts_in_plane(:,m)/np1
        enddo
! connect points in the clockwise order
!  two points are in clockwise order when:
!   1) scalar product of vectors is maximum AND
!   2) vector product is collinear to the normal of the plane 
        do n=1,np1
! first vector
          q1(:)=pts_in_plane(:,n)-q3(:)
          a1=-100.d0
          do m=1,np1
            if (m.ne.n) then
! second vector
              q2(:)=pts_in_plane(:,m)-q3(:)
! scalar product              
              a3=q1(1)*q2(1)+q1(2)*q2(2)+q1(3)*q2(3)
! vector product             
              q4(1)=q1(2)*q2(3)-q2(2)*q1(3)
              q4(2)=-(q1(1)*q2(3)-q2(1)*q1(3)) 
              q4(3)=q1(1)*q2(2)-q2(1)*q1(2)
! scalar product with the norm to check collinearity              
              a2=q4(1)*norm(1)+q4(2)*norm(2)+q4(3)*norm(3)
              if (a2.gt.1d-4.and.a3.gt.a1) then
                a1=a3
                j=m
              endif
            endif
          enddo !m
          conn(1,n)=pts_idx(n)
          conn(2,n)=pts_idx(j)
        enddo !n
! find point with minmal index from which to start building the face
!  this in needed for easy comparison of polygons
        j1=100
        do n=1,np1
          if (conn(1,n).lt.j1) j1=conn(1,n)
        enddo
! build the face of the polyhedra
        face(0)=np1
        n=1
        do while (n.le.np1)
          do m=1,np1
            if (conn(1,m).eq.j1) then
              face(n)=j1
              j1=conn(2,m)
              n=n+1
            endif
          enddo
        enddo
! add new connections to the list
        do n=1,np1
          l1=.true.
          do m=1,nc
            if ((connections(1,m).eq.conn(1,n).and.connections(2,m).eq.conn(2,n)).or. &
                (connections(1,m).eq.conn(2,n).and.connections(2,m).eq.conn(1,n)) ) then
              l1=.false.
            endif
          enddo !m
          if (l1) then
            nc=nc+1
            connections(1,nc)=conn(1,n)
            connections(2,nc)=conn(2,n)
          endif
        enddo !n
! add new face
        l1=.true.
        do m=1,nf
          if ((faces(0,m).eq.face(0)).and. &
              (all(faces(1:faces(0,m),m).eq.face(1:face(0))))) then
            l1=.false.
          endif
        enddo !m
        if (l1) then
          nf=nf+1
          faces(:,nf)=face(:)
        endif
      endif
    enddo
  enddo
enddo
return
end


