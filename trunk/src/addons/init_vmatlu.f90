subroutine init_vmatlu
use modmain
use modldapu
implicit none
integer ispn,m,l,is,ias,ia
logical file_is_found
real(8), allocatable :: orbital_occupancy(:,:)

inquire(file="dmat.in",exist=file_is_found)
if (.not.file_is_found) return
open(50,file="dmat.in",form="FORMATTED",status="OLD")

dmatlu=zzero
do is=1,nspecies
  l=llu(is)
  if (l.gt.0) then
    allocate(orbital_occupancy(2*l+1,nspinor))
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      read(50,*)((orbital_occupancy(m,ispn),m=1,2*l+1),ispn=1,nspinor)
      do ispn=1,nspinor
        do m=1,2*l+1
          dmatlu(l**2+m,l**2+m,ispn,ispn,ias)=dcmplx(orbital_occupancy(m,ispn),0.d0)
        enddo
      enddo
    enddo
    deallocate(orbital_occupancy)
  endif
enddo 
do ias=1,natmtot
  do ispn=1,nspinor
    call unimtrxt(lmmaxlu,yrlm_lps(1,1,ias),dmatlu(1,1,ispn,ispn,ias))
  enddo
enddo
call genvmatlu
if (wproc) call writeldapu
return
end subroutine
