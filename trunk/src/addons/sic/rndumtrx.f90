subroutine rndumtrx(n,mtrx)
implicit none
integer, intent(in) :: n
complex(8), intent(out) :: mtrx(n,n)
!
integer i,j
real rn(2)
real(8), allocatable :: eval(:)
complex(8), allocatable :: zm(:,:)
!
allocate(eval(n))
allocate(zm(n,n))
do i=1,n
  do j=1,n
    call random_number(rn)
    zm(i,j)=dcmplx(rn(1),rn(2))
  enddo
enddo
do i=1,n
  do j=1,n
    mtrx(i,j)=zm(i,j)+dconjg(zm(j,i))
  enddo
enddo
call diagzhe(n,mtrx,eval)
deallocate(eval,zm)
return
end subroutine
