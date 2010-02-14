subroutine init_g_q_gq(ivq0m,lmaxexp,lmmaxexp)
use modmain
implicit none
integer, intent(in) :: ivq0m(3)
integer, intent(in) :: lmaxexp
integer, intent(in) :: lmmaxexp
integer vgq0l(3)
integer i,j,ig

! q-vector in lattice coordinates
do i=1,3
  vq0l(i)=1.d0*(ivq0m(i))/ngridk(i)+1d-12
enddo
! find G-vector which brings q0 to first BZ
vgq0l(:)=floor(vq0l(:))
! reduce q0 vector to first BZ
vq0rl(:)=vq0l(:)-vgq0l(:)
! get Cartesian coordinates of q-vector and reduced q-vector
call r3mv(bvec,vq0l,vq0c)
call r3mv(bvec,vq0rl,vq0rc)

call getngvecme(vgq0l)

if (allocated(lr_vgq0c)) deallocate(lr_vgq0c)
allocate(lr_vgq0c(3,ngvecme))
if (allocated(lr_gq0)) deallocate(lr_gq0)
allocate(lr_gq0(ngvecme))
if (allocated(lr_tpgq0)) deallocate(lr_tpgq0)
allocate(lr_tpgq0(2,ngvecme))
if (allocated(lr_sfacgq0)) deallocate(lr_sfacgq0)
allocate(lr_sfacgq0(ngvecme,natmtot))
if (allocated(lr_ylmgq0)) deallocate(lr_ylmgq0)
allocate(lr_ylmgq0(lmmaxexp,ngvecme))

! check if we have enough G-shells to bring q-vector back to first BZ
do ig=1,ngvecme
  if (sum(abs(vgq0l(:)-ivg(:,ig+gvecme1-1))).eq.0) then
    lr_igq0=ig+gvecme1-1
    goto 20
  endif
enddo
write(*,*)
write(*,'("Error(init_g_q_gq): no G-vector to reduce q-vector to first BZ")')
write(*,'("  ngvecme : ",I4)')ngvecme
write(*,'("  vq0l : ",3G18.10)')vq0l
write(*,'("  vgq0l : ",3G18.10)')vgq0l
write(*,*)
call pstop
20 continue

! generate G+q' vectors, where q' is reduced q-vector
i=0
do ig=gvecme1,gvecme2
  i=i+1
  lr_vgq0c(:,i)=vgc(:,ig)+vq0rc(:)
! get spherical coordinates and length of G+q'
  call sphcrd(lr_vgq0c(:,i),lr_gq0(i),lr_tpgq0(:,i))
! generate spherical harmonics for G+q'
  call genylm(lmaxexp,lr_tpgq0(:,i),lr_ylmgq0(:,i))
enddo

! generate structure factor for G+q' vectors
call gensfacgp(ngvecme,lr_vgq0c,ngvecme,lr_sfacgq0)

return
end