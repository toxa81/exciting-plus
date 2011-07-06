subroutine getnufr
use modmain
implicit none

integer ltmp(0:lolmax)
integer is,ilo,l,nlomaxl,io
!
! TODO: fix for the case when some species have APW and some LAPW
!
nufr=apwordmax
! find maximum number of local orbitals over all l-channels
nlomaxl=0
do is=1,nspecies
  ltmp=0
  do ilo=1,nlorb(is)
    ltmp(lorbl(ilo,is))=ltmp(lorbl(ilo,is))+1
    nufr(lorbl(ilo,is),is)=nufr(lorbl(ilo,is),is)+1
  enddo
  do l=0,lolmax
    if (l.le.lmaxapw) then
      nlomaxl=max(nlomaxl,ltmp(l))
    endif
  enddo
enddo
nufrmax=apwordmax+nlomaxl
nlufr=0
do is=1,nspecies
  do l=0,lmaxapw
    do io=1,nufr(l,is)
      nlufr(is)=nlufr(is)+1
    enddo
  enddo
enddo
nlufrmax=maxval(nlufr)
return
end
