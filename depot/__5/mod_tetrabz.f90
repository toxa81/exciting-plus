module mod_tetrabz

real(8) qpts(3,26)
integer nbzpt
real(8) bzpts(3,200)
integer ntetra
real(8), allocatable :: tetras(:,:,:)

contains

logical function find_bz_vertex(iq1,iq2,iq3,v)
implicit none
integer, intent(in) :: iq1
integer, intent(in) :: iq2
integer, intent(in) :: iq3
real(8), intent(out) :: v(3)
!
integer i
real(8) a(3,3),ai(3,3),b(3),d0,d1
real(8), external :: r3mdet
!
find_bz_vertex=.false.
a(1,:)=qpts(:,iq1)/2.d0
a(2,:)=qpts(:,iq2)/2.d0
a(3,:)=qpts(:,iq3)/2.d0

if (abs(r3mdet(a)).lt.1d-10) return

do i=1,3
  b(i)=dot_product(a(i,:),a(i,:))
enddo
call r3minv(a,ai)
call r3mv(ai,b,v)

d0=dot_product(v,v)
do i=1,26
  b(:)=v(:)-qpts(:,i)
  d1=dot_product(b,b)
  if (d0.gt.(d1+1d-10)) return
enddo
find_bz_vertex=.true.
return
end function

subroutine tetrabz_split
use modmain
implicit none
real(8), allocatable :: tetras_new(:,:,:)
allocate(tetras_new(3,4,ntetra*8))


return
end subroutine


subroutine tetrabz_fill
use modmain
implicit none
integer i1,i2,i3,i,j,n1,n2,j1,j2,nfp,nfacets,itet
integer facet_(20),facet(20),facets(0:20,40)
real(8) v(3),q0(3),q1(3),q2(3),q3(3),norm(3),prod,max_prod,vol
logical tfound
j=0
do i1=-1,1
  do i2=-1,1
    do i3=-1,1
      if (.not.(i1.eq.0.and.i2.eq.0.and.i3.eq.0)) then
        j=j+1
        qpts(:,j)=i1*bvec(:,1)+i2*bvec(:,2)+i3*bvec(:,3)
      endif
    enddo
  enddo
enddo

nbzpt=0
do i1=1,24
  do i2=i1+1,25
    do i3=i2+1,26
      if (find_bz_vertex(i1,i2,i3,v)) then
        tfound=.false.
        do j=1,nbzpt
          if (sum(abs(v(:)-bzpts(:,j))).lt.1d-10) tfound=.true.
        enddo
        if (.not.tfound) then
          nbzpt=nbzpt+1
          bzpts(:,nbzpt)=v(:)
        endif
      endif
    enddo
  enddo
enddo

nfacets=0
do i1=1,nbzpt-2
  do i2=i1+1,nbzpt-1
    do i3=i2+1,nbzpt
      q1(:)=bzpts(:,i1)-bzpts(:,i3)
      q2(:)=bzpts(:,i2)-bzpts(:,i3)
      call r3cross(q1,q2,norm)
      q0=0.d0
      j1=nbzpt+1
      n1=0
      n2=0
      nfp=0
      do j=1,nbzpt
        q3(:)=bzpts(:,j)-bzpts(:,i3)
        prod=dot_product(q3,norm)
        if (abs(prod).lt.1d-10) then
          nfp=nfp+1
          facet_(nfp)=j
          if (j.lt.j1) j1=j
          q0(:)=q0(:)+bzpts(:,j)
        else if (prod.ge.1d-10) then
          n1=n1+1
        else 
          n2=n2+1
        endif
      enddo
      if (n1.eq.0.or.n2.eq.0) then
        if (n2.eq.0)  norm=-1.d0*norm
        q0=q0/nfp
        do n1=1,nfp
          q1(:)=bzpts(:,j1)-q0(:)
          max_prod=-100.d0
          do n2=1,nfp
            q2(:)=bzpts(:,facet_(n2))-q0(:)
            call r3cross(q1,q2,q3)
            prod=dot_product(q1,q2)
            if (prod.gt.max_prod.and.dot_product(q3,norm).gt.0.d0) then
              max_prod=prod
              j2=facet_(n2)
            endif
          enddo
          facet(n1)=j1
          j1=j2
        enddo !n1
        tfound=.false.
        do j=1,nfacets
          if (facets(0,j).eq.nfp) then
            if (all(facets(1:nfp,j).eq.facet(1:nfp))) tfound=.true.
          endif
        enddo
        if (.not.tfound) then
          nfacets=nfacets+1
          facets(0,nfacets)=nfp
          facets(1:nfp,nfacets)=facet(1:nfp)
        endif
      endif
    enddo !i3
  enddo !i2
enddo !i1
write(*,*)"number of facets : ",nfacets
ntetra=0
do j=1,nfacets
  nfp=facets(0,j)
  !do i=1,nfp 
  !  write(*,*)facets(i,j),facets(mod(i,nfp)+1,j)
  !enddo
  !write(*,*)
  ntetra=ntetra+nfp
enddo
allocate(tetras(3,4,ntetra))
itet=0
do j=1,nfacets
  nfp=facets(0,j)
  q0=0.d0
  do i=1,nfp
    q0=q0+bzpts(:,facets(i,j))/nfp
  enddo
  do i=1,nfp
    i1=facets(i,j)
    i2=facets(mod(i,nfp)+1,j)
    itet=itet+1
    tetras(:,1,itet)=0.d0
    tetras(:,2,itet)=q0(:)
    tetras(:,3,itet)=bzpts(:,i1)
    tetras(:,4,itet)=bzpts(:,i2)
  enddo
enddo
vol=0.d0
do itet=1,ntetra
  q1(:)=tetras(:,1,itet)-tetras(:,4,itet)
  q2(:)=tetras(:,2,itet)-tetras(:,4,itet)
  q3(:)=tetras(:,3,itet)-tetras(:,4,itet)
  call r3cross(q1,q2,q0)
  vol=vol+abs(dot_product(q0,q3))/6.d0
enddo
write(*,*)"volume=",vol
 
    







stop
     
      


return
end subroutine

end module
