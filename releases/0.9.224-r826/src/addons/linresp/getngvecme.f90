subroutine getngvecme(vgq0l)
use modmain
implicit none
integer, intent(in) :: vgq0l(3)
integer ngsh,gshq0,i,j
integer, allocatable :: igishell(:)
integer, allocatable :: ishellng(:,:)

! find G-shell for a given q-vector
allocate(igishell(ngvec))
allocate(ishellng(ngvec,2))
call getgshells(ngsh,igishell,ishellng)
if (vgq0l(1).lt.intgv(1,1).or.vgq0l(1).gt.intgv(1,2).or.&
    vgq0l(2).lt.intgv(2,1).or.vgq0l(2).gt.intgv(2,2).or.&
    vgq0l(3).lt.intgv(3,1).or.vgq0l(3).gt.intgv(3,2)) then
  write(*,*)
  write(*,'("Error(getngvecme): Gq-vector is outside of boundaries")')
  write(*,'("                   (q-vector is too large)")')
  write(*,'("  Gq : ",3I8)')vgq0l
  write(*,'("  boundaries : ",2I5,",",2I5,",",2I5)')intgv(1,:),intgv(2,:),&
    intgv(3,:)
  call pstop
endif
    
gshq0=igishell(ivgig(vgq0l(1),vgq0l(2),vgq0l(3)))
if (wproc) then
  write(150,*)
  write(150,'("G-shell of a given q-vector : ",I4)')gshq0
endif
if (gshq0.lt.gshme1) then
  if (wproc) then
    write(150,*)
    write(150,'("Warning: minimum number of G-shells was changed from ",&
      &I4," to ",I4)')gshme1,gshq0
  endif
  gshme1=gshq0
endif
if (gshq0.gt.gshme2) then
  if (wproc) then
    write(150,*)
    write(150,'("Warning: maximum number of G-shells was changed from ",&
      &I4," to ",I4)')gshme2,gshq0
  endif
  gshme2=gshq0
endif
! test if G-shells are closed
i=ishellng(gshme1,2)
j=ishellng(gshme2,2)
if (abs(gc(i)-gc(i+1)).lt.epslat.or.abs(gc(j)-gc(j+1)).lt.epslat) then
  write(*,*)
  write(*,'("Error(getngvecme): G-shells are not closed")')
  write(*,*)
  call pstop
endif
if (gshme1.eq.1) then
  gvecme1=1
else
  gvecme1=ishellng(gshme1-1,2)+1
endif
gvecme2=ishellng(gshme2,2)
ngvecme=gvecme2-gvecme1+1
if (scalar_chi) then
  if (wproc) then
    write(150,*)
    write(150,'("Scalar calculation")')
  endif
  gvecme1=ivgig(vgq0l(1),vgq0l(2),vgq0l(3))
  gvecme2=gvecme1
!  gshme1=-1
!  gshme2=gshme1
  ngvecme=1
endif
deallocate(igishell)
deallocate(ishellng)
return
end