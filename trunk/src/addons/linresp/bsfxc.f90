subroutine bsfxc(iq,chi0m,fxckrnl)
use modmain
use mod_expigqr
use mod_linresp
use mod_addons_q
implicit none
integer, intent(in) :: iq
complex(8), intent(in) :: chi0m(ngq(iq),ngq(iq))
complex(8), intent(inout) :: fxckrnl(ngq(iq),ngq(iq))
integer iter,ig,ig1,ig2
real(8) tdiff
complex(8), allocatable :: krnl(:,:)
complex(8), allocatable :: eps(:,:)
complex(8), allocatable :: chi(:,:)
complex(8), allocatable :: fxckrnl1(:,:)

allocate(krnl(ngq(iq),ngq(iq)))
allocate(eps(ngq(iq),ngq(iq)))
allocate(chi(ngq(iq),ngq(iq)))
allocate(fxckrnl1(ngq(iq),ngq(iq)))

do iter=1,200
! restore the full kernel
  krnl=fxckrnl
  do ig=1,ngq(iq)
    krnl(ig,ig)=krnl(ig,ig)+vhgq(ig,iq)
  enddo
! compute epsilon=1-chi0*(v+fxc) 
  eps=zzero
  do ig=1,ngq(iq)
    eps(ig,ig)=zone
  enddo
  call zgemm('N','N',ngq(iq),ngq(iq),ngq(iq),-zone,chi0m,ngq(iq),krnl,&
    ngq(iq),zone,eps,ngq(iq))
! invert epsilon matrix
  call invzge(eps,ngq(iq))
! compute chi=epsilon^-1 * chi0
  call zgemm('N','N',ngq(iq),ngq(iq),ngq(iq),zone,eps,ngq(iq),chi0m,&
    ngq(iq),zzero,chi,ngq(iq))
! compute screend Coulomb potential vscr=vbare+vbare*chi*vbare
  fxckrnl1=zzero
  do ig=1,ngq(iq)
    fxckrnl1(ig,ig)=vhgq(ig,iq)
  enddo
  do ig1=1,ngq(iq)
    do ig2=1,ngq(iq)
      fxckrnl1(ig1,ig2)=fxckrnl1(ig1,ig2)+vhgq(ig1,iq)*chi(ig1,ig2)*vhgq(ig2,iq)
    enddo
  enddo


  krnl=zzero
  do ig=1,ngq(iq)
    krnl(ig,ig)=krnl(ig,ig)+vhgq(ig,iq)
  enddo
  call zgemm('N','N',ngq(iq),ngq(iq),ngq(iq),zone,chi0m,ngq(iq),krnl,&
    ngq(iq),zzero,eps,ngq(iq))
! invert epsilon matrix
  call invzge(eps,ngq(iq))
! compute chi=epsilon^-1 * w
  call zgemm('N','N',ngq(iq),ngq(iq),ngq(iq),zone,fxckrnl1,ngq(iq),eps,&
    ngq(iq),zzero,krnl,ngq(iq))
    fxckrnl1=krnl

! scale kernel by eps0
   !do ig1=1,ngq(iq)
   ! do ig2=1,ngq(iq)
   !   fxckrnl1(ig1,ig2)=fxckrnl1(ig1,ig2)/(chi0m(iig0q,iig0q)*vhgq(iig0q,iq))
   !   !if (ig1.ne.ig2) fxckrnl(ig1,ig2)=zzero
   ! enddo
  !enddo
  tdiff=0.d0
  do ig1=1,ngq(iq)
    do ig2=1,ngq(iq)
      tdiff=tdiff+abs(fxckrnl(ig1,ig2)-fxckrnl1(ig1,ig2))
    enddo
  enddo
  fxckrnl=0.75d0*fxckrnl+0.25d0*fxckrnl1
  if (tdiff.lt.1d-8) goto 10
enddo
write(*,'("Warning(bsfx): total difference : ",G18.10)')tdiff
10 continue
deallocate(krnl,eps,chi,fxckrnl1)
return
end subroutine
