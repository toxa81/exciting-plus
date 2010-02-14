subroutine sumchi0(ikloc,w,chi0w)
use modmain
implicit none
integer, intent(in) :: ikloc
complex(8), intent(in) :: w
complex(8), intent(out) :: chi0w(ngvecme,ngvecme)

logical l1
complex(8) zt1
integer i,ist1,ist2
integer, parameter :: bs=128
integer, parameter :: chi0summation=4
integer nb,sz1,offs,ik,jk
integer ib1,ib2,j1,j2,ig
logical, allocatable :: l2(:)
complex(8), allocatable :: wt(:)
logical, external :: bndint

ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
jk=idxkq(1,ik)
offs=nmegqblhloc(2,ikloc)
allocate(l2(nmegqblhloc(1,ikloc)))
allocate(wt(nmegqblhlocmax))
l2=.false.
do i=1,nmegqblhloc(1,ikloc)
  ist1=bmegqblh(1,i+offs,ikloc)
  ist2=bmegqblh(2,i+offs,ikloc)
! default : include all interband transitions         
  l1=.true.
! for cRPA : don't include bands in energy window [crpa_e1,crpa_e2]
  if (crpa) then
    if (bndint(ist1,lr_evalsvnr(ist1,ik),crpa_e1,crpa_e2).and. &
        bndint(ist2,lr_evalsvnr(ist2,jk),crpa_e1,crpa_e2)) l1=.false.
  endif
  if (l1) then
    if (abs(lr_occsvnr(ist1,ik)-lr_occsvnr(ist2,jk)).gt.1d-5) then
      wt(i)=(lr_occsvnr(ist1,ik)-lr_occsvnr(ist2,jk))/(lr_evalsvnr(ist1,ik) - &
        lr_evalsvnr(ist2,jk)+w)
      l2(i)=.true.
    endif
  else
    wt(i)=zzero
  endif
enddo !i

if (chi0summation.eq.1) then
  do i=1,nmegqblhloc(1,ikloc)
    if (l2(i)) then
      call zgerc(ngvecme,ngvecme,wt(i),megqblh(i,:,ikloc),1,megqblh(i,:,ikloc),1, &
        chi0w(1,1),ngvecme)
    endif
  enddo
endif
if (chi0summation.eq.2) then
! number of blocks
  nb=ngvecme/bs
! remaining size
  sz1=mod(ngvecme,bs)
  ib1=1
  do j1=1,nb
    ib2=1
    do j2=1,nb
      do i=1,2 !i1,i2
        if (l2(i)) then
          call zgerc(bs,bs,wt(i),megqblh(ib1,i,ikloc),1,megqblh(ib2,i,ikloc),1, &
            chi0w(ib1,ib2),ngvecme)
        endif
      enddo !i
      ib2=ib2+bs
    enddo !j2
    ib1=ib1+bs
  enddo !j1
! remaining part
  if (sz1.ne.0) then
    ib1=1
    do j1=1,nb
      do i=1,2 !i1,i2
        if (l2(i)) then
          call zgerc(bs,sz1,wt(i),megqblh(ib1,i,ikloc),1,megqblh(nb*bs+1,i,ikloc),1, &
            chi0w(ib1,nb*bs+1),ngvecme)
          call zgerc(sz1,bs,wt(i),megqblh(nb*bs+1,i,ikloc),1,megqblh(ib1,i,ikloc),1, &
            chi0w(nb*bs+1,ib1),ngvecme)
        endif
      enddo !i
      ib1=ib1+bs
    enddo !j1
    do i=1,2 !i1,i2
      if (l2(i)) then
        call zgerc(sz1,sz1,wt(i),megqblh(nb*bs+1,i,ikloc),1,megqblh(nb*bs+1,i,ikloc),1, &
          chi0w(nb*bs+1,nb*bs+1),ngvecme)
      endif
    enddo !i
  endif
endif
if (chi0summation.eq.4) then
  do ig=1,ngvecme
    megqblh2(:,ig)=dconjg(megqblh(:,ig,ikloc))*wt(:)
  enddo
  call zgemm('T','N',ngvecme,ngvecme,nmegqblhloc(1,ikloc),zone,&
    megqblh(1,1,ikloc),nmegqblhlocmax,megqblh2(1,1),nmegqblhlocmax,&
    zone,chi0w(1,1),ngvecme)
endif

deallocate(l2)
deallocate(wt)
return
end
