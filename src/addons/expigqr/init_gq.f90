subroutine init_gq(iq,lmaxexp,lmmaxexp,tg0q)
use modmain
use mod_addons_q
implicit none
integer, intent(in) :: iq
integer, intent(in) :: lmaxexp
integer, intent(in) :: lmmaxexp
logical, intent(in) :: tg0q
integer vgq0l(3)
real(8) v1(3)
integer i,ig
real(8) t1

ngvecme=ngq(iq)
!if (wproc) then
!  write(*,'(" iq : ",I6,5X,"ngvecme : ",I6)')iq,ngvecme
!endif
if (allocated(tpgq)) deallocate(tpgq)
allocate(tpgq(2,ngvecme))
if (allocated(sfacgq)) deallocate(sfacgq)
allocate(sfacgq(ngvecme,natmtot))
if (allocated(ylmgq)) deallocate(ylmgq)
allocate(ylmgq(lmmaxexp,ngvecme))
! check if we have enough G-shells to bring q-vector back to first BZ
if (tg0q) then
  do ig=1,ngvecme
    if (igqig(ig,iq).eq.ig0q(iq)) then
      iig0q=ig
      goto 20
    endif
  enddo
  write(*,*)
  write(*,'("Error(init_gq): no G-vector to reduce q-vector to first BZ")')
  v1=vqc(:,iq)+vgc(:,ig0q(iq))
  write(*,'(" q (mesh coord.) : ",3I6)')vqm(:,iq)
  write(*,'(" G0 : ",3I6)')ivg(:,ig0q(iq))
  write(*,'(" |G0+q| : ",G18.10)')sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
  write(*,'(" gqmax : ",G18.10)')gqmax
  write(*,'(" gqsh : ",I4)')gqsh
  write(*,'(" ngvecme: ",I4)')ngvecme    
  write(*,*)
  call pstop
  20 continue
endif

! generate spherical harmonics for G+q vectors
do ig=1,ngvecme
! get spherical coordinates and length of G+q
  call sphcrd(vgqc(:,ig,iq),t1,tpgq(:,ig))
! generate spherical harmonics for G+q'
  call genylm(lmaxexp,tpgq(:,ig),ylmgq(:,ig))
enddo

! generate structure factor for G+q vectors
call gensfacgp(ngvecme,vgqc(1,1,iq),ngvecme,sfacgq)

return
end