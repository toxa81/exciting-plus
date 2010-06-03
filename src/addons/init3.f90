subroutine init3
use modmain
implicit none
integer ia,is,lm,l,m
if (allocated(rylm)) deallocate(rylm)
allocate(rylm(16,16))
if (allocated(yrlm)) deallocate(yrlm)
allocate(yrlm(16,16))
if (allocated(rylm_lps)) deallocate(rylm_lps)
allocate(rylm_lps(16,16,natmtot))
if (allocated(yrlm_lps)) deallocate(yrlm_lps)
allocate(yrlm_lps(16,16,natmtot))
call genshmat
if (allocated(ias2is)) deallocate(ias2is)
allocate(ias2is(natmtot))
if (allocated(ias2ia)) deallocate(ias2ia)
allocate(ias2ia(natmtot))
do is=1,nspecies
  do ia=1,natoms(is)
    ias2is(idxas(ia,is))=is
    ias2ia(idxas(ia,is))=ia
  end do
end do
call getatmcls
if (allocated(nufr)) deallocate(nufr)
allocate(nufr(0:lmaxvr,nspecies))
call getnufr(lmaxvr)
if (allocated(ufr)) deallocate(ufr)
allocate(ufr(nrmtmax,0:lmaxvr,nufrmax,natmcls))
if (allocated(ufrp)) deallocate(ufrp)
allocate(ufrp(0:lmaxvr,nufrmax,nufrmax,natmcls))
if (allocated(lm2l)) deallocate(lm2l)
allocate(lm2l(lmmaxapw))
if (allocated(lm2m)) deallocate(lm2m)
allocate(lm2m(lmmaxapw))
do l=0,lmaxapw
  do m=-l,l
    lm2l(idxlm(l,m))=l
    lm2m(idxlm(l,m))=m
  end do
end do
if (wannier) call wann_init
return
end