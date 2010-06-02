subroutine getatmcls
use modmain
implicit none
integer ias,is,ia1,ia2,i
natmcls=0
if (allocated(ias2ic)) deallocate(ias2ic)
allocate(ias2ic(natmtot))
ias2ic=0
do is=1,nspecies
  do ia1=1,natoms(is)
    if (ias2ic(idxas(ia1,is)).eq.0) then
      natmcls=natmcls+1
      do ia2=1,natoms(is) 
        if (eqatoms(ia1,ia2,is)) then
          ias2ic(idxas(ia2,is))=natmcls
        endif
      enddo
    endif
  enddo
enddo
if (allocated(ic2ias)) deallocate(ic2ias)
allocate(ic2ias(natmcls))
do i=1,natmcls
  do ias=1,natmtot
    if (ias2ic(ias).eq.i) then
      ic2ias(i)=ias 
      exit
    endif
  enddo
enddo
return
end
