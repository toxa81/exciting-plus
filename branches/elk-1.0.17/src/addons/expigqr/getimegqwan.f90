subroutine getimegqwan(lall)
use modmain
implicit none
logical, intent(in) :: lall
integer n,n1,i,ias,jas,j,j1
logical l1
logical, external :: wann_diel

if (nwann_include.eq.0) then
  nwann_include=nwann
  if (allocated(iwann_include)) deallocate(iwann_include)
  allocate(iwann_include(nwann))
  do j=1,nwann
    iwann_include(j)=j
  enddo
endif

call getnghbr(megqwan_mindist,megqwan_maxdist)
! get maximum possible number of WF transitions
nmegqwanmax=0
do n=1,nwann
  ias=iwann(1,n)
  do i=1,nnghbr(ias)
    do n1=1,nwann
      jas=iwann(1,n1)
      if (jas.eq.inghbr(1,i,ias)) then
        nmegqwanmax=nmegqwanmax+nwannias(jas)
      endif
    enddo
  enddo
enddo
if (allocated(imegqwan)) deallocate(imegqwan)
allocate(imegqwan(5,nmegqwanmax))
imegqwan=0
nmegqwan=0   
do j=1,nwann_include
  n=iwann_include(j)
  ias=iwann(1,n)
  do i=1,nnghbr(ias)
    do j1=1,nwann_include
      n1=iwann_include(j1)
      jas=iwann(1,n1)
      if (jas.eq.inghbr(1,i,ias)) then
        l1=.false.
! for integer occupancy numbers take only transitions between occupied and empty bands
        if (wann_diel().and.(abs(wann_occ(n)-wann_occ(n1)).gt.1d-8)) l1=.true.
! for fractional occupancies or other cases take all transitions
        if (.not.wann_diel().or.lall) l1=.true.
        if (l1) then
          nmegqwan=nmegqwan+1
          imegqwan(1,nmegqwan)=n
          imegqwan(2,nmegqwan)=n1
          imegqwan(3:5,nmegqwan)=inghbr(3:5,i,ias)
        endif 
      endif
    enddo
  enddo
enddo
megqwan_tlim(1,1)=minval(imegqwan(3,:))
megqwan_tlim(2,1)=maxval(imegqwan(3,:))
megqwan_tlim(1,2)=minval(imegqwan(4,:))
megqwan_tlim(2,2)=maxval(imegqwan(4,:))
megqwan_tlim(1,3)=minval(imegqwan(5,:))
megqwan_tlim(2,3)=maxval(imegqwan(5,:))
if (allocated(idxmegqwan)) deallocate(idxmegqwan)
allocate(idxmegqwan(nwann,nwann,megqwan_tlim(1,1):megqwan_tlim(2,1),&
  megqwan_tlim(1,2):megqwan_tlim(2,2),megqwan_tlim(1,3):megqwan_tlim(2,3)))
idxmegqwan=-100
do i=1,nmegqwan
  idxmegqwan(imegqwan(1,i),imegqwan(2,i),imegqwan(3,i),imegqwan(4,i),&
    imegqwan(5,i))=i
enddo
if (allocated(itmegqwan)) deallocate(itmegqwan)
allocate(itmegqwan(3,nmegqwan))
ntmegqwan=0
do i=1,nmegqwan
  l1=.true.
  do j=1,ntmegqwan
    if (all(itmegqwan(:,j).eq.imegqwan(3:5,i))) l1=.false.
  enddo
  if (l1) then
    ntmegqwan=ntmegqwan+1
    itmegqwan(:,ntmegqwan)=imegqwan(3:5,i)
  endif
enddo
return
end