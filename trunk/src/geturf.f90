subroutine geturf
use modmain
implicit none
! arguments

integer is,ia,ias,l,io,ilo
integer ordl(0:lmaxvr)

urf=0.d0
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    ordl=0
! apw functions
    do l=0,lmaxvr
      do io=1,apword(l,is)
        ordl(l)=ordl(l)+1
        urf(1:nrmt(is),l,ordl(l),ias)=apwfr(1:nrmt(is),1,io,l,ias)
      enddo
    enddo
! lo functions
    do ilo=1,nlorb(is)
      l=lorbl(ilo,is)
      if (l.le.lmaxvr) then
        ordl(l)=ordl(l)+1
        urf(1:nrmt(is),l,ordl(l),ias)=lofr(1:nrmt(is),1,ilo,ias)
      endif
    enddo
  enddo
enddo

return
end

