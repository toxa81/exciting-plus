subroutine genrlm_lcs
use modmain
implicit none
integer lmax,lmmax
integer isym,lspl,lm,lm1,lm2,lm3,l,n,m1,m2
integer is,ia,ias,i,j
complex(8), allocatable :: zm1(:,:)
complex(8), allocatable :: zm2(:,:)
complex(8), allocatable :: ulat(:,:,:)
complex(8), allocatable :: a(:,:),b(:,:)
complex(8), allocatable :: h(:,:)
real(8), allocatable :: mtrx(:,:)
real(8), allocatable :: eval(:)
complex(8), allocatable :: mtrx1(:,:)


call init0

lmax=3
lmmax=(lmax+1)**2

allocate(ulat(lmmax,lmmax,nsymlat))
allocate(a(lmmax,lmmax),b(lmmax,lmmax))
allocate(h(lmmax,lmmax))
allocate(zm1(lmmax,lmmax))
allocate(zm2(lmmax,lmmax))





! construct (l,m) rotation matrix for each lattice symmetry
a(:,:)=0.d0
do i=1,lmmax
  a(i,i)=1.d0
end do
do isym=1,nsymlat
  call rotzflm(symlatc(:,:,isym),lmax,lmmax,lmmax,a,ulat(:,:,isym))
end do
! set up quasi-random symmetric matrix H
h(:,:)=0.d0
do l=0,lmax
  n=2*l+1
  lm=idxlm(l,-l)
  do i=lm,lm+n-1
    do j=i,lm+n-1
      h(i,j)=dble((i+0.1d0)*(j+0.2d0))
      h(j,i)=h(i,j)
    end do
  end do
end do
zm1=zzero
zm2=zzero
do lm1=1,lmmax
  do lm2=1,lmmax
    do lm3=1,lmmax
      zm2(lm1,lm2)=zm2(lm1,lm2)+dconjg(yrlm(lm1,lm3))*h(lm3,lm2)
    enddo
  enddo
enddo
do lm1=1,lmmax
  do lm2=1,lmmax
    do lm3=1,lmmax
      zm1(lm1,lm2)=zm1(lm1,lm2)+zm2(lm1,lm3)*yrlm(lm2,lm3)
    enddo
  enddo
enddo
h=zm1


! loop over species and atoms
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! symmetrise H with site symmetries
    b(:,:)=0.d0
    do isym=1,nsymsite(ias)
! spatial rotation element in lattice point group
      lspl=lsplsyms(isym,ias)
! apply lattice symmetry as U*H*conjg(U')
      call zgemm('N','N',lmmax,lmmax,lmmax,zone,ulat(:,:,lspl),lmmax,h,lmmax, &
       zzero,a,lmmax)
      call zgemm('N','C',lmmax,lmmax,lmmax,zone,a,lmmax,ulat(:,:,lspl),lmmax, &
       zone,b,lmmax)
    end do !isym

    zm1=zzero
    zm2=zzero
    do lm1=1,lmmax
      do lm2=1,lmmax
        do lm3=1,lmmax
          zm2(lm1,lm2)=zm2(lm1,lm2)+dconjg(rylm(lm1,lm3))*b(lm3,lm2)
        enddo
      enddo
    enddo
    do lm1=1,lmmax
      do lm2=1,lmmax
        do lm3=1,lmmax
          zm1(lm1,lm2)=zm1(lm1,lm2)+zm2(lm1,lm3)*rylm(lm2,lm3)
        enddo
      enddo
    enddo

    write(*,*)'ias=',ias
    l=2
    allocate(mtrx(2*l+1,2*l+1))
    allocate(eval(2*l+1))
    do m1=-l,l
      do m2=-l,l
        mtrx(m1+l+1,m2+l+1)=dreal(zm1(idxlm(l,m1),idxlm(l,m2)))
      enddo
    enddo
    call diagdsy(2*l+1,mtrx,eval)
    write(*,*)'LCS matrix'
    do m1=1,2*l+1
      write(*,'(7G18.10)')(mtrx(m1,m2),m2=1,2*l+1)
    enddo
    write(*,*)'Transpose of matrix'
    do m1=1,2*l+1
      write(*,'(7G18.10)')(mtrx(m2,m1),m2=1,2*l+1)
    enddo
    write(*,'("  eval : ",7G18.10)')eval(:)
    deallocate(mtrx,eval)
  end do
end do
deallocate(ulat,a,b,h,zm1,zm2)
return
end