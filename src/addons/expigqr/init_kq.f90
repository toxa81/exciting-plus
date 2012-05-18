subroutine init_kq(iq)
use modmain
use mod_addons_q
use mod_wannier
use mod_expigqr
implicit none
integer, intent(in) :: iq
integer ik,jk
integer ivg1(3)
real(8) vkql(3)
real(8), external :: r3taxi

! find k+q and reduce them to first BZ (this is required to utilize the 
!   periodical property of Bloch-states: |k>=|k+K>, where K is any vector 
!   of the reciprocal lattice)
if (allocated(idxkq)) deallocate(idxkq)
allocate(idxkq(3,nkptnr))
do ik=1,nkptnr
! k+q vector
  vkql(:)=vklnr(:,ik)+vql(:,iq)+1d-12
! K vector
  ivg1(:)=floor(vkql(:))
! reduced k+q vector: k'=k+q-K
  vkql(:)=vkql(:)-ivg1(:)
! search for index of reduced k+q vector 
  do jk=1,nkptnr
    if (r3taxi(vklnr(:,jk),vkql).lt.epslat) then
      idxkq(1,ik)=jk
      goto 10
    endif
  enddo
  write(*,*)
  write(*,'("Error(init_kq): index of reduced k+q point is not found")')
  write(*,'(" index of k-point: ",I4)')ik
  write(*,'(" K-vector: ",3I4)')ivg1
  write(*,'(" reduced k+q vector: ",3G18.10)')vkql
  write(*,'(" check original q-vector coordinates")')
  write(*,*)
  call pstop
10 continue
  idxkq(2,ik)=ivgig(ivg1(1),ivg1(2),ivg1(3))
! find k-q index
  vkql(:)=vklnr(:,ik)-vql(:,iq)+1d-12
! K vector
  ivg1(:)=floor(vkql(:))
! reduced k-q vector: k'=k-q-K
  vkql(:)=vkql(:)-ivg1(:)
! search for index of reduced k-q vector 
  do jk=1,nkptnr
    if (r3taxi(vklnr(:,jk),vkql).lt.epslat) then
      idxkq(3,ik)=jk
      goto 20
    endif
  enddo
  write(*,*)
  write(*,'("Error(init_kq): index of reduced k-q point is not found")')
  write(*,'(" index of k-point: ",I4)')ik
  write(*,'(" K-vector: ",3I4)')ivg1
  write(*,'(" reduced k-q vector: ",3G18.10)')vkql
  write(*,'(" check original q-vector coordinates")')
  write(*,*)
  call pstop
20 continue 
enddo
! check index correctness
do ik=1,nkptnr
! k+q point
  jk=idxkq(1,ik)
! check index of k-q point at jk 
  if (idxkq(3,jk).ne.ik) then
    write(*,*)
    write(*,'("Error(init_kq): k+q and k-q indices are not consitent")')
    write(*,'(" index of k-point : ",I4)')ik
    write(*,'(" index of k+q point : ",I4)')jk
    write(*,'(" index of k-q point at jk : ",I4)')idxkq(3,jk)
    call pstop
  endif
enddo
return
end
