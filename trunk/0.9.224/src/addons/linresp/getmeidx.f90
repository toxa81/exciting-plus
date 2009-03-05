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
! possible candidate for charge response
        if (ispn1.eq.ispn2.and.lrtype.eq.0) then
          if ((ispn1.eq.spin_me.or.spin_me.eq.3).and.ldocc) laddme=.true.
        endif
! for magnetic response
        if (ispn1.ne.ispn2.and.lrtype.eq.1) then
          if ((ispn1.eq.spin_me.or.spin_me.eq.3).and.ldocc) laddme=.true.
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

