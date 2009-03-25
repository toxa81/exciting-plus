subroutine getmeidx(req,occsvnr)
use modmain
implicit none
! arguments
logical, intent(in) :: req
real(8), intent(in) :: occsvnr(nstsv,nkptnr)

integer band1,band2
integer i,ik,jk,istfv1,istfv2,ispn1,ispn2,ist1,ist2,ikloc
logical laddme,ldocc
real(8) d1
logical l11,l12,l21,l22
integer, external :: iknrglob2

if (bndme1.eq.-1) then
  band1=1
  band2=nstfv
else
  band1=bndme1
  band2=bndme2
endif
if (req.and.wproc) then
  write(150,*)
  write(150,'("Band interval: ",2I4)')band1,band2
endif
if (req) nmemax=0
do ikloc=1,nkptnr_loc
  ik=iknrglob2(ikloc,mpi_x(1))
  jk=idxkq(1,ik)
  i=0
  do ispn1=1,nspinor
    do ispn2=1,nspinor
      do istfv1=band1,band2
      do istfv2=band1,band2
        ist1=istfv1+(ispn1-1)*nstfv
        ist2=istfv2+(ispn2-1)*nstfv
        d1=occsvnr(ist1,ik)-occsvnr(ist2,jk)
        ldocc=abs(d1).gt.1d-10
        laddme=.false.
	if (ldocc) then
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
            ime(3,i,ikloc)=ispn1
            docc(i,ikloc)=d1
          endif
        endif
      enddo !istfv2
      enddo !istfv1
    enddo !ispn2
  enddo !ispn1
  if (.not.req) nme(ikloc)=i
  if (req) nmemax=max(nmemax,i)
enddo !ikloc

return
end

