subroutine getnghbr(mindist,maxdist)
use modmain
implicit none
real(8), intent(in) :: mindist
real(8), intent(in) :: maxdist

integer i1,i2,i3,ias,ia,is,jas,ja,js,i,n,j
real(8) a(3),v1(3),d1,v2(3),v3(3),dst
integer llim(3),tmp(5)
real(8) vrc(3),vrl(3)
integer ntr(3)

! find lattice limits
llim=0
do i=1,3
  vrc=0.d0
  vrc(i)=maxdist
  call getntr(vrc,ntr,vrl)
  do j=1,3
    llim(j)=max(abs(ntr(j))+1,llim(j))
  enddo
enddo

if (allocated(nnghbr)) deallocate(nnghbr)
allocate(nnghbr(natmtot))
nnghbr=0
if (allocated(inghbr)) deallocate(inghbr)
allocate(inghbr(5,natmtot*(2*llim(1)+1)*(2*llim(2)+1)*(2*llim(3)+1),natmtot))
inghbr=0

do ias=1,natmtot
  ia=ias2ia(ias)
  is=ias2is(ias)
  do jas=1,natmtot
    ja=ias2ia(jas)
    js=ias2is(jas)
    do i1=-llim(1),llim(1)
      do i2=-llim(2),llim(2)
        do i3=-llim(3),llim(3)
          v1(:)=atposc(:,ja,js)+i1*avec(:,1)+i2*avec(:,2)+i3*avec(:,3)-&
            atposc(:,ia,is)
          d1=sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
          if (d1.ge.mindist.and.d1.le.maxdist) then
            nnghbr(ias)=nnghbr(ias)+1
            n=nnghbr(ias)
            inghbr(1,n,ias)=jas
            inghbr(2,n,ias)=int(d1*1000000)
            inghbr(3,n,ias)=i1
            inghbr(4,n,ias)=i2
            inghbr(5,n,ias)=i3
          endif
        enddo !i3
      enddo !i2
    enddo !i1
  enddo !jas
! sort by distance    
  do i1=1,nnghbr(ias)-1
    do i2=i1+1,nnghbr(ias)
      if (inghbr(2,i1,ias).gt.inghbr(2,i2,ias)) then
        tmp(:)=inghbr(:,i1,ias)
        inghbr(:,i1,ias)=inghbr(:,i2,ias)
        inghbr(:,i2,ias)=tmp(:)
      endif
    enddo
  enddo
enddo    
return
end