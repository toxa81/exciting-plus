subroutine getnfr(lmax)
use modmain
implicit none
integer, intent(in) :: lmax

integer ltmp(0:lolmax)
integer is,ilo,l,nlomaxl
nfr=apwordmax
! find maximum number of local orbitals over all l-channels
nlomaxl=0
do is=1,nspecies
  ltmp=0
  do ilo=1,nlorb(is)
    ltmp(lorbl(ilo,is))=ltmp(lorbl(ilo,is))+1
    nfr(lorbl(ilo,is),is)=nfr(lorbl(ilo,is),is)+1
  enddo
  do l=0,lolmax
    if (l.le.lmax) then
      nlomaxl=max(nlomaxl,ltmp(l))
    endif
  enddo
enddo
nfrmax=apwordmax+nlomaxl
return
end
