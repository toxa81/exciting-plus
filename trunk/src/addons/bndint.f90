logical function bndint(j,e,e1,e2)
use modmain
use mod_wannier
implicit none
integer, intent(in) :: j
real(8), intent(in) :: e
real(8), intent(in) :: e1
real(8), intent(in) :: e2

logical leint1,leint2,l1,l2
integer n1,n2

leint1=.true.
leint2=.true.
if (e1.gt.0.d0) then
  n1=int(e1)
  if (abs(dble(n1)-e1).lt.1d-5) leint1=.false.
endif
if (e2.gt.0.d0) then
  n2=int(e2)
  if (abs(dble(n2)-e2).lt.1d-5) leint2=.false.
endif

if (leint1) then
  l1=(e.ge.e1)
else
 if ((ndmag.eq.3).or.ldisentangle) then
    l1=(j.ge.n1)
  else
    l1=((mod(j-1,nstfv)+1).ge.n1)
  endif
endif
if (leint2) then
  l2=(e.le.e2)
else
 if ((ndmag.eq.3).or.ldisentangle) then
    l2=(j.le.n2)
  else
    l2=((mod(j-1,nstfv)+1).le.n2)
  endif
endif
bndint=(l1.and.l2)
!leint=.true.
!if (abs(n1-e1).lt.1d-5.and.abs(n2-e2).lt.1d-5) then
!  leint=.false.
!endif
!if (leint) then
!  bndint=(e.ge.e1.and.e.le.e2)
!else
!  if (ncmag.or.ldisentangle) then
!    bndint=(j.ge.n1.and.j.le.n2)
!  else
!    bndint=((mod(j-1,nstfv)+1).ge.n1.and.(mod(j-1,nstfv)+1).le.n2)
!  endif
!endif
return
end function

