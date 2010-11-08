subroutine svdchi0(chi0wan,iq,iw,w)
use modmain
use mod_addons_q

! arguments
complex(8), intent(in) :: chi0wan(nmegqwan,nmegqwan)
integer, intent(in) :: iq
integer, intent(in) :: iw
real(8), intent(in) :: w
! local variables
real(8), allocatable :: A(:,:)
real(8), allocatable :: S(:)
real(8), allocatable :: U(:,:)
real(8), allocatable :: VT(:,:)
real(8), allocatable :: tmp1(:)
real(8), allocatable :: tmp2(:)

integer i,j,iwf1,iwf2
integer lda,ldu,ldvt,lwork,info,lwmax
real(8), allocatable :: work(:)

integer nlow
real(8) t1

character*100 fname
!character*10 c1

integer ias,jas,lm1,lm2,n,n1,vtl(3)
real(8) vtc(3),vrc(3)
character*20 c1,c2,c3,c4
character, parameter :: orb(4)=(/'s','p','d','f'/)
character*200 str

! do allocations
lwmax = 100000
allocate(A(nmegqwan,nmegqwan))
allocate(S(nmegqwan))
allocate(U(nmegqwan,nmegqwan))
allocate(VT(nmegqwan,nmegqwan))
allocate(work(lwmax))
allocate(tmp1(nmegqwan))
allocate(tmp2(nmegqwan))

! make copy of chi0wan
A = abs(chi0wan)

! lapack nonsense
lwork = lwmax
lda = nmegqwan
ldu = nmegqwan
ldvt = nmegqwan

!write(c1,'(I6)')iw
!fname="matrixA__"//trim(adjustl(c1))//".txt"
!open(210,file=trim(fname),form='formatted',status='replace')
!if (wproc) then
!  write(210,'(F15.5)')w
!  do i=1,nmegqwan
!    write(210,'(100F15.5)')A(i,:)
!  enddo
!endif
!close(210)


!if (wproc) then
!  print *, ''
!  print *, ''
!  print *, ''
!  print *, 'A matrix:'
!  print *, '------------------------------------------------------------'
!  do i=1,nmegqwan
!    print '(100F15.5)',A(i,:)
!  enddo
!endif

call dgesvd('A','A',nmegqwan,nmegqwan,A(1,1),lda,S(1),U(1,1),ldu,VT(1,1),ldvt,work,&
  lwork,info)

nlow = -1
do i=1,nmegqwan
  t1 = S(1)/S(i)
  if (t1.lt.4.0) then
    nlow = i
  endif
enddo

!if (wproc) then
!  print *, ''
!  print *, ''
!  print *, ''
!  print *, 'S matrix:'
!  print *, '------------------------------------------------------------'
!  do i=1,nmegqwan
!    t1 = S(1)/S(i)
!    if (t1.lt.4.0) then
!      nlow = i
!    endif
!    write(*,'(I6," singular value:  ",E15.5)')i,S(i)
!  enddo
!
!  write(*,'("The largest singular value: ", E15.5)'), S(1)
!  write( *,'("omega: ", F15.5)'), w
!  write( *,'("nlow: ", I6)'), nlow
!  
!  
!  print *, ''
!  print *, ''
!  print *, 'U matrix:'
!  print *, '------------------------------------------------------------'
!  do i=1,nlow
!    print '(100E15.5)',U(:,i)
!  enddo
!  print *, ''
!  print *, ''
!  print *, ''
!endif

! loop over particle-hole pairs
do i=1,nmegqwan
  ! loop over singular values
  do j=1,nmegqwan
    tmp1(i)=tmp1(i)+S(j)*abs(U(i,j))
  enddo
enddo

if (wproc) then
  write(c1,'(I6)')iw
  fname="sing2_values__"//trim(adjustl(c1))//".txt"
  open(210,file=trim(fname),form='formatted',status='replace')
  write(210,'(F15.5)')w
  do i=1,nmegqwan
    n=imegqwan(1,i)
    n1=imegqwan(2,i)
    write(210,'(I4," ---> ",I4,"     : ",F15.5)')n,n1,tmp1(i)
  enddo
  close(210)
endif

if (wproc) then
  write(c1,'(I6)')iw
  fname="sing_values__"//trim(adjustl(c1))//".txt"
  open(210,file=trim(fname),form='formatted',status='replace')
  write(210,'(F15.5)')w
  do i=1,nlow
      str="           "
      write(210,'(A)')trim(str)
      write(210,'(A)')trim(str)
      write(210,'(A)')trim(str)
      write(210,'(A)')trim(str)
      write(c1,'(G18.10)')S(i)
      write(210,'(I4)')i
      str="           SV="//trim(adjustl(c1))
      write(210,'(A)')
    do j=1,nmegqwan
      if (U(j,i).gt.0.1) then
        n=imegqwan(1,j)
        n1=imegqwan(2,j)
        vtl=imegqwan(3:5,j)
        ias=iwann(1,n)
        lm1=iwann(2,n)
        jas=iwann(1,n1)
        lm2=iwann(2,n1)
        vtc(:)=avec(:,1)*vtl(1)+avec(:,2)*vtl(2)+avec(:,3)*vtl(3)
        vrc(:)=vtc(:)+atposc(:,ias2ia(jas),ias2is(jas))-&
                      atposc(:,ias2ia(ias),ias2is(ias))       
        write(c1,'(I6)')ias2ia(ias)
        write(c2,'(I1)')lm2m(lm1)+lm2l(lm1)+1
        c3="("//trim(spsymb(ias2is(ias)))//trim(adjustl(c1))//"-"//&
          orb(lm2l(lm1)+1)//trim(adjustl(c2))//")"
        write(c1,'(I6)')ias2ia(jas)
        write(c2,'(I1)')lm2m(lm2)+lm2l(lm2)+1
        c4="("//trim(spsymb(ias2is(jas)))//trim(adjustl(c1))//"-"//&
          orb(lm2l(lm2)+1)//trim(adjustl(c2))//")"
        write(c1,'("j : ",I4)')j
        str=trim(adjustl(c1))
        write(c1,'(I6)')n
        str=trim(str)//"   n="//trim(adjustl(c1))//" "//trim(adjustl(c3))
        write(c1,'(I6)')n1
        str=trim(str)//" ->  n'="//trim(adjustl(c1))//" "//trim(adjustl(c4))
        write(c1,'(I6)')vtl(1)
        write(c2,'(I6)')vtl(2)
        write(c3,'(I6)')vtl(3)
        str=trim(str)//" T=("//trim(adjustl(c1))//" "//trim(adjustl(c2))//" "//&
          trim(adjustl(c3))//")"
        write(c1,'(F12.4)')vrc(1)
        write(c2,'(F12.4)')vrc(2)
        write(c3,'(F12.4)')vrc(3)
        str=trim(str)//" R=("//trim(adjustl(c1))//" "//trim(adjustl(c2))//" "//&
          trim(adjustl(c3))//")"
        write(c1,'(F12.4)')sqrt(sum(vrc(:)**2))
        str=trim(str)//" D="//trim(adjustl(c1))
        write(c1,'(G18.10)')U(j,i)
        write(210,'(A)')trim(str)
        str="           coeff="//trim(adjustl(c1))
        write(210,'(A)')trim(str)
      endif
    enddo
  enddo
  close(210)
endif

! deallocate
deallocate(A)
deallocate(S)
deallocate(U)
deallocate(VT)
deallocate(tmp1)
deallocate(tmp2)

end

