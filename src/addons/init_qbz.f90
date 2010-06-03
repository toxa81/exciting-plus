subroutine init_qbz(lgamma)
use modmain
implicit none
logical, intent(in) :: lgamma
integer j,i1,i2,i3

if (allocated(ivq0m_list)) deallocate(ivq0m_list)
nvq0=nkptnr-1
j=0  
if (lgamma) then
  nvq0=nvq0+8
  j=8
endif
allocate(ivq0m_list(3,nvq0))
ivq0m_list=0
do i1=0,ngridk(1)-1
  do i2=0,ngridk(2)-1
    do i3=0,ngridk(3)-1
      if (.not.(i1.eq.0.and.i2.eq.0.and.i3.eq.0)) then
        j=j+1
        ivq0m_list(:,j)=(/i1,i2,i3/)
      endif
    enddo
  enddo
enddo
if (lgamma) call init_q0gamma
return
end
