subroutine bsfxc(iq,chi0m,fxckrnl)
use modmain
use mod_expigqr
use mod_linresp
use mod_addons_q
implicit none
integer, intent(in) :: iq
complex(8), intent(in) :: chi0m(ngvecme,ngvecme)
complex(8), intent(inout) :: fxckrnl(ngvecme,ngvecme)
integer iter,ig,ig1,ig2
real(8) tdiff
complex(8), allocatable :: krnl(:,:)
complex(8), allocatable :: eps(:,:)
complex(8), allocatable :: chi(:,:)

allocate(krnl(ngvecme,ngvecme))
allocate(eps(ngvecme,ngvecme))
allocate(chi(ngvecme,ngvecme))

do iter=1,50
! restore the full kernel
  krnl=fxckrnl
  do ig=1,ngvecme
    krnl(ig,ig)=krnl(ig,ig)+vhgq(ig,iq)
  enddo
! compute epsilon=1-chi0*(v+fxc) 
  eps=zzero
  do ig=1,ngvecme
    eps(ig,ig)=zone
  enddo
  call zgemm('N','N',ngvecme,ngvecme,ngvecme,-zone,chi0m,ngvecme,krnl,&
    ngvecme,zone,eps,ngvecme)
! invert epsilon matrix
  call invzge(eps,ngvecme)
! compute chi=epsilon^-1 * chi0
  call zgemm('N','N',ngvecme,ngvecme,ngvecme,zone,eps,ngvecme,chi0m,&
    ngvecme,zzero,chi,ngvecme)
! compute screend Coulomb potential vscr=vbare+vbare*chi*vbare
  fxckrnl=zzero
  do ig=1,ngvecme
    fxckrnl(ig,ig)=vhgq(ig,iq)
  enddo
  do ig1=1,ngvecme
    do ig2=1,ngvecme
      fxckrnl(ig1,ig2)=fxckrnl(ig1,ig2)+vhgq(ig1,iq)*chi(ig1,ig2)*vhgq(ig2,iq)
    enddo
  enddo
! scale kernel by eps0
   do ig1=1,ngvecme
    do ig2=1,ngvecme
      fxckrnl(ig1,ig2)=fxckrnl(ig1,ig2)/(chi0m(iig0q,iig0q)*vhgq(iig0q,iq))
      if (ig1.ne.ig2) fxckrnl(ig1,ig2)=zzero
    enddo
  enddo
  tdiff=0.d0
  do ig1=1,ngvecme
    do ig2=1,ngvecme
      if (ig1.eq.ig2) then
        tdiff=tdiff+abs(fxckrnl(ig1,ig1)-krnl(ig1,ig1)+vhgq(ig1,iq))
      else
        tdiff=tdiff+abs(fxckrnl(ig1,ig2)-krnl(ig1,ig2))
      endif
    enddo
  enddo
  if (tdiff.lt.1d-8) goto 10
enddo
write(*,'("Warning(bsfx): total difference : ",G18.10)')tdiff
10 continue
deallocate(krnl,eps,chi)
return
end subroutine
