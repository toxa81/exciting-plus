subroutine getmeidx(req,occsvnr,evalsvnr)
use modmain
implicit none
! arguments
logical, intent(in) :: req
real(8), intent(in) :: occsvnr(nstsv,nkptnr)
real(8), intent(in) :: evalsvnr(nstsv,nkptnr)

integer band1,band2
integer i,ik,jk,ist1,ist2,ikloc
logical laddme,ldocc
real(8) d1
logical l11,l12,l21,l22,le1,le2
integer, external :: iknrglob2

if (req.and.wproc) then
  write(150,*)
  write(150,'("Band interval (Ha) : ",2F8.2)')lr_e1,lr_e2
endif
if (req) nmemax=0
do ikloc=1,nkptnr_loc
  ik=iknrglob2(ikloc,mpi_x(1))
  jk=idxkq(1,ik)
  i=0
  do ist1=1,nstsv
    do ist2=1,nstsv
      le1=evalsvnr(ist1,ik).ge.lr_e1.and.evalsvnr(ist1,ik).le.lr_e2
      le2=evalsvnr(ist2,jk).ge.lr_e1.and.evalsvnr(ist2,jk).le.lr_e2
      d1=occsvnr(ist1,ik)-occsvnr(ist2,jk)
      ldocc=abs(d1).gt.1d-10
      laddme=.false.
      if (ldocc.and.le1.and.le2) then
        if (.not.spinpol) then
	  laddme=.true.
	else
	  l11=spinor_ud(1,ist1,ik).eq.1.and.spinor_ud(1,ist2,jk).eq.1
	  l12=spinor_ud(1,ist1,ik).eq.1.and.spinor_ud(2,ist2,jk).eq.1
	  l21=spinor_ud(2,ist1,ik).eq.1.and.spinor_ud(1,ist2,jk).eq.1
	  l22=spinor_ud(2,ist1,ik).eq.1.and.spinor_ud(2,ist2,jk).eq.1
	  if (lrtype.eq.0.and.(l11.or.l22)) laddme=.true.
	  if (lrtype.eq.1.and.(l12.or.l21)) laddme=.true.
	endif
      endif
      if (laddme) then
        i=i+1
        if (.not.req) then
          ime(1,i,ikloc)=ist1
          ime(2,i,ikloc)=ist2
          ime(3,i,ikloc)=1 
          docc(i,ikloc)=d1
        endif
      endif
    enddo
  enddo
  if (.not.req) nme(ikloc)=i
  if (req) nmemax=max(nmemax,i)
enddo !ikloc

return
end

