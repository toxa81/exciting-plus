subroutine getufr
use modmain
implicit none
integer is,ias,l,io,ilo,ic
integer ordl(0:lmaxapw)
ufr=0.d0
do ic=1,natmcls
  ias=ic2ias(ic)
  is=ias2is(ias)
! apw functions
  do l=0,lmaxapw
    do io=1,apword(l,is)
      ufr(1:nrmt(is),l,io,ic)=apwfr(1:nrmt(is),1,io,l,ias)
    enddo
    ordl(l)=apword(l,is)
  enddo
! lo functions
  do ilo=1,nlorb(is)
    l=lorbl(ilo,is)
    if (l.le.lmaxapw) then
      ordl(l)=ordl(l)+1
      ufr(1:nrmt(is),l,ordl(l),ic)=lofr(1:nrmt(is),1,ilo,ias)
    endif
  enddo
enddo
return
end

