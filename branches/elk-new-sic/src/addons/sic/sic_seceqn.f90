subroutine sic_seceqn(ikloc,evecfv,evecsv)
use modmain
use mod_lf
use mod_nrkp
implicit none
! arguments
integer, intent(in) :: ikloc
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(out) :: evecsv(nstsv,nstsv)
! local variables
integer ispn,ist,n,j
complex(8), allocatable :: wann_ufv(:,:,:)

allocate(wann_ufv(nwann,nstfv,nspinor))
wann_ufv=zzero
do ispn=1,nspinor
  do ist=1,nstfv
    do n=1,nwann
      do j=1,nstsv
        wann_ufv(n,ist,ispn)=wann_ufv(n,ist,ispn)+&
          wann_c(n,j,ikloc)*evecsv(ist+(ispn-1)*nstfv,j)
      enddo
    enddo
  enddo
enddo
call sic_hunif(ikloc,wann_ufv,evecfv,evecsv)
deallocate(wann_ufv)
return
end