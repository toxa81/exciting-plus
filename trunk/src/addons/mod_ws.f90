module mod_ws

real(8) ws_vec(3,3)
real(8) ws_vecinv(3,3)

! nearest neigboring points
real(8) ws_nnpts(3,26)
integer ws_nvertex
real(8) ws_vertex(3,512)

contains

logical function ws_find_vertex(i1,i2,i3,v)
implicit none
integer, intent(in) :: i1
integer, intent(in) :: i2
integer, intent(in) :: i3
integer i
real(8), intent(in) :: v(3)
real(8) a(3,3),ai(3,3),b(3),d0,d1
real(8), external :: r3mdet
!
ws_find_vertex=.false.
a(1,:)=ws_nnpts(:,i1)/2.d0
a(2,:)=ws_nnpts(:,i2)/2.d0
a(3,:)=ws_nnpts(:,i3)/2.d0

if (abs(r3mdet(a)).lt.1d-10) return

call r3minv(a,ai)
do i=1,3
  b(i)=dot_product(a(i,:),a(i,:))
enddo
call r3mv(ai,b,v)

d0=dot_product(v,v)
do i=1,26
  b(:)=ws_nnpts(:,i)-v(:)
  d1=dot_product(b,b)
  if (d0.gt.(d1+1d-10)) return
enddo
ws_find_vertex=.true.
return
end function

subroutine ws_reduce(v,rmax,t)
implicit none
real(8), intent(inout) :: v(3)
real(8), optional, intent(in) :: rmax
real(8), optional, intent(out) :: t(3)
integer i,j
real(8) d,d1,v1(3)
! V_reduced = V_orig + T
j=0
d=dot_product(v,v)
do i=1,26
  v1(:)=v(:)+ws_nnpts(:,i)
  d1=dot_product(v1,v1)
  if (d1.lt.d) then
    d=d1
    j=i
  endif
enddo
if (present(t)) t=0
if (j.ne.0) then
  v=v+ws_nnpts(:,j)
  if (present(t)) t(:)=ws_nnpts(:,j)
endif
if (present(rmax)) then
  if (sum(v(:)**2).gt.(rmax**2+1d-8)) then
    write(*,'("Error(ws_reduce): outside of rmax")')
    write(*,'("  rmax : ",G18.10)')rmax
    write(*,'("   x   : ",3G18.10)')v(:)
    write(*,'("  |x|  : ",G18.10)')sqrt(sum(v(:)**2))
    call pstop
  endif
endif
return
end subroutine

subroutine ws_init(a)
implicit none
real(8), intent(in) :: a(3,3)
real(8) v(3)
integer j,i1,i2,i3
logical tfound
!
ws_vec=a
call r3minv(ws_vec,ws_vecinv)
j=0
do i1=-1,1
  do i2=-1,1
    do i3=-1,1
      if (.not.(i1.eq.0.and.i2.eq.0.and.i3.eq.0)) then
        j=j+1
        ws_nnpts(:,j)=i1*a(:,1)+i2*a(:,2)+i3*a(:,3)
      endif
    enddo
  enddo
enddo

ws_nvertex=0
do i1=1,24
  do i2=i1+1,25
    do i3=i2+1,26
      if (ws_find_vertex(i1,i2,i3,v)) then
        tfound=.false.
        do j=1,ws_nvertex
          if (sum(abs(v(:)-ws_vertex(:,j))).lt.1d-10) tfound=.true.
        enddo
        if (.not.tfound) then
          ws_nvertex=ws_nvertex+1
          ws_vertex(:,ws_nvertex)=v(:)
        endif
      endif
    enddo
  enddo
enddo
      
      
      

return
end subroutine

logical function ws_vector_inside(v)
implicit none
real(8), intent(in) :: v(3)
!
integer i
real(8) d,d1,v1(3)
d=dot_product(v,v)
ws_vector_inside=.true.
do i=1,26
  v1(:)=v(:)+ws_nnpts(:,i)
  d1=dot_product(v1,v1)
  if (d1.lt.d) then
    ws_vector_inside=.false.
    return
  endif
enddo
return
end function

end module
