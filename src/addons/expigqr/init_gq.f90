subroutine init_gq(iq,lmaxexp,lmmaxexp,tg0q)
use modmain
use mod_addons_q
use mod_wannier
use mod_expigqr
implicit none
integer, intent(in) :: iq
integer, intent(in) :: lmaxexp
integer, intent(in) :: lmmaxexp
logical, intent(in) :: tg0q
real(8) v1(3)
integer ig
real(8) t1

if (allocated(tpgq)) deallocate(tpgq)
allocate(tpgq(2,ngq(iq)))
if (allocated(sfacgq)) deallocate(sfacgq)
allocate(sfacgq(ngq(iq),natmtot))
if (allocated(ylmgq)) deallocate(ylmgq)
allocate(ylmgq(lmmaxexp,ngq(iq)))
! check if we have enough G-shells to bring q-vector back to first BZ
iig0q=-1
if (tg0q) then
  do ig=1,ngq(iq)
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
  write(*,'(" ngq(iq): ",I4)')ngq(iq)    
  write(*,*)
  call pstop
  20 continue
endif

! generate spherical harmonics for G+q vectors
do ig=1,ngq(iq)
! get spherical coordinates and length of G+q
  call sphcrd(vgqc(:,ig,iq),t1,tpgq(:,ig))
! generate spherical harmonics for G+q'
  call genylm(lmaxexp,tpgq(:,ig),ylmgq(:,ig))
enddo

! generate structure factor for G+q vectors
call gensfacgp(ngq(iq),vgqc(1,1,iq),ngq(iq),sfacgq)

return
end