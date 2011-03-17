logical function bndint(j,e,e1,e2)
use modmain
use mod_wannier
implicit none
integer, intent(in) :: j
real(8), intent(in) :: e
real(8), intent(in) :: e1
real(8), intent(in) :: e2

logical leint
integer n1,n2

leint=.true.
n1=int(e1)
n2=int(e2)
if (abs(n1-e1).lt.1d-5.and.abs(n2-e2).lt.1d-5) then
  leint=.false.
endif
if (leint) then
  bndint=(e.ge.e1.and.e.le.e2)
else
  if (ncmag.or.ldisentangle) then
    bndint=(j.ge.n1.and.j.le.n2)
  else
    bndint=((mod(j-1,nstfv)+1).ge.n1.and.(mod(j-1,nstfv)+1).le.n2)
  endif
endif
return
end
