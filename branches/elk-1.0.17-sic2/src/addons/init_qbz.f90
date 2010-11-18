subroutine init_qbz(lgamma,nvq0_)
use modmain
use mod_addons_q
implicit none
logical, intent(in) :: lgamma
integer, intent(in) :: nvq0_
integer j,i1,i2,i3
real(8) t1

if (.not.(nvq0_.eq.1.or.nvq0_.eq.8)) then
  write(*,*)
  write(*,'("Error(init_qbz): number of q0 vectors must be 1 or 8")')
  call pstop
endif
if (allocated(vqm)) deallocate(vqm)
nvq=nkptnr-1
nvq0=nvq0_
j=0  
if (lgamma) then
  nvq=nvq+nvq0
  j=nvq0
endif
allocate(vqm(3,nvq))
vqm=0
do i1=0,ngridk(1)-1
  do i2=0,ngridk(2)-1
    do i3=0,ngridk(3)-1
      if (.not.(i1.eq.0.and.i2.eq.0.and.i3.eq.0)) then
        j=j+1
        vqm(:,j)=(/i1,i2,i3/)
      endif
    enddo
  enddo
enddo
if (lgamma) call init_q0
! for unscreened U 
if (lgamma.and.nvq0.eq.1) then
  t1=0.d0
  do j=1,8
    t1=t1+aq0(j)/dot_product(vq0c(:,j),vq0c(:,j))
  enddo
  aq0(1)=t1
endif
return
end
