subroutine writebz
use modmain
implicit none
real(8) pt1(3,26),p0(3),pt2(3,100),pt3(3,200)
integer pl(200)
integer i1,i2,i3,i,j,m,npt,nc,nf,n1,nptot
integer connections(2,100)
integer faces(0:10,100)
logical vbz,l1
real(8) b(3)

call init0
call init1

! connect center point with neigbours
m=1
do i1=-1,1
  do i2=-1,1
    do i3=-1,1
      if (.not.(i1.eq.0.and.i2.eq.0.and.i3.eq.0)) then
        pt1(:,m)=i1*bvec(:,1)+i2*bvec(:,2)+i3*bvec(:,3)
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

nptot=npt+nc+nf
pt3(:,1:npt)=pt2(:,1:npt)
do i=1,nc
  pt3(:,npt+i)=0.5d0*(pt2(:,connections(1,i))+pt2(:,connections(2,i)))
enddo
do i=1,nf
  p0=0.d0
  do m=1,faces(0,i)
    p0(:)=p0(:)+pt2(:,faces(m,i))
  enddo
  pt3(:,npt+nc+i)=p0/faces(0,i)
enddo

! count points starting from 0
faces(1:10,:)=faces(1:10,:)-1
connections(:,:)=connections(:,:)-1

n1=0
do i=1,nf
  n1=n1+faces(0,i)
enddo

open(150,file='bz.dx',status='replace',form='formatted')
! labels of the points
write(150,'("object 1 class array type string rank 1 shape 4 items ",&
  &I3," data follows")')nptot
pl=1
do m=1,nptot
  do j=1,nptot
! skip inverse points
    if (sum(abs(pt3(:,m)+pt3(:,j))).lt.1d-8) pl(j)=0
  enddo
! gat lattice coordinates of the point
  call r3mv(binv,pt3(:,m),b)
  if (pl(m).eq.1) then
    write(150,'("""",I3,""""," # cart=",3F12.6," lat=",3F12.6)')m,pt3(:,m),b
  else
    write(150,'("""   """," # cart=",3F12.6," lat=",3F12.6)')pt3(:,m),b    
  endif  
enddo
write(150,'("attribute ""dep"" string ""positions""")')
! color of points
write(150,'("object 2 class array type float rank 1 shape 3 items ",&
  &I3," data follows")')nptot
do i=1,nptot
  write(150,'(3F12.6)')0.d0,0.d0,1.d0
enddo
write(150,'("attribute ""dep"" string ""positions""")')
! positions
write(150,'("object 3 class array type float rank 1 shape 3 items ",&
  &I3," data follows")')nptot
do i=1,nptot
  write(150,'(3F12.6)')pt3(:,i)
enddo
write(150,'("attribute ""dep"" string ""positions""")')
! connections
write(150,'("object 4 class array type int rank 1 shape 2 items ",&
  &I3," data follows")')nc
do i=1,nc
  write(150,'(2I4)')connections(:,i)
enddo
write(150,'("attribute ""element type"" string ""lines""")')
write(150,'("attribute ""ref"" string ""positions""")')
! faces
write(150,'("object 5 class array type int rank 0 items ",&
  &I3," data follows")')n1
do i=1,nf
  write(150,'(10I4)')faces(1:faces(0,i),i)
enddo
write(150,'("attribute ""ref"" string ""positions""")')
write(150,'("object 6 class array type int rank 0 items ",&
  &I3," data follows")')nf
i=0
j=1
do while (i.lt.n1)
  write(150,'(I3)')i
  i=i+faces(0,j)
  j=j+1
enddo
write(150,'("attribute ""ref"" string ""edges""")')
write(150,'("object 7 class array type int rank 0 items ",&
  &I3," data follows")')nf
do i=1,nf
  write(150,'(I3)')i-1
enddo
write(150,'("attribute ""ref"" string ""loops""")')
write(150,'("object ""bz"" class field")')
write(150,'("component ""data"" value 1")')
write(150,'("component ""colors"" value 2")')
write(150,'("component ""positions"" value 3")')
write(150,'("component ""connections"" value 4")')
write(150,'("component ""edges"" value 5")') 
write(150,'("component ""loops"" value 6")')
write(150,'("component ""faces"" value 7")')
write(150,'("end")')
close(150)

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
real(8) a(3,3),b(3),a1(3,3)
real(8) det0
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
if (abs(det0).lt.1d-10) then
  p0(:)=10000.d0
  return
endif

! b-vector for the system of linear equations
b(1)=(q1(1)**2+q1(2)**2+q1(3)**2)
b(2)=(q2(1)**2+q2(2)**2+q2(3)**2)
b(3)=(q3(1)**2+q3(2)**2+q3(3)**2)

! invert 3x3 matrix
call r3minv(a,a1)
! x=a^{-1}*b
call r3mv(a1,b,p0)

return
end


subroutine is_vector_in_bz(p0,pt1,f)
implicit none
real(8), intent(in) :: p0(3)
real(8), intent(in) :: pt1(3,26)
logical, intent(out) :: f

real(8) d0,d1
integer i

! squared distance from the center of BZ to the point p0
d0=p0(1)**2+p0(2)**2+p0(3)**2
f=.true.

do i=1,26
! squared distance from the nearest neigbour lattice point to the point p0
  d1=(p0(1)-pt1(1,i))**2+(p0(2)-pt1(2,i))**2+(p0(3)-pt1(3,i))**2
! check the main property of BZ: for any point inside BZ distance to center 
!   is smaller than distance to neighbouring lattice points
  if (d0.gt.(d1+1d-10)) f=.false.
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
      if (all(fpt(:).gt.-1d-10).or.all(fpt(:).lt.1d-10)) then
! norm should point out of BZ (to the side of the plane without points)
        if (all(fpt(:).gt.-1d-10)) norm(:)=-norm(:)
! find all points lying in the plane; this are the points, satisfying the plane
!  equation n*(x-p)=0
        np1=0
        do m=1,np
          if (abs(fpt(m)).lt.1d-10) then
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
          j=-100
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
              if (a2.gt.1d-10.and.a3.gt.a1) then
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


