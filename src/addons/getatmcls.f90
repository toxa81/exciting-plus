subroutine getatmcls
use modmain
implicit none
integer ias,is,ia1,ia2,i,ic
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
allocate(natoms_in_class(natmcls))
natoms_in_class=0
do ias=1,natmtot
  natoms_in_class(ias2ic(ias))=natoms_in_class(ias2ic(ias))+1
enddo
if (iproc.eq.0) then
  open(200,file="CLASS.OUT",form="FORMATTED",status="REPLACE")
  do ias=1,natmtot
    write(200,'("ias : ",I4,"   ic : ",I4)')ias,ias2ic(ias)
  enddo
  write(200,*)
  do ic=1,natmcls
    write(200,'("ic : ",I4,"   ias : ",I4," natoms : ",I4)')ic,ic2ias(ic),natoms_in_class(ic)  
  enddo
  close(200)
endif
return
end
