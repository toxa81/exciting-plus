subroutine getmeidx(req,occsvnr,evalsvnr)
use modmain
#ifdef _MPI_
use mpi
#endif
implicit none
! arguments
logical, intent(in) :: req
real(8), intent(in) :: occsvnr(nstsv,nkptnr)
real(8), intent(in) :: evalsvnr(nstsv,nkptnr)

integer i,ik,jk,ist1,ist2,ikloc,n
logical laddme,ldocc
real(8) d1,min_e12,min_e1_wann,max_e2_wann
logical l11,l12,l21,l22,le1,le2,lwann,le1w,le2w
integer, external :: iknrglob2
logical, external :: bndint
logical, external :: wann_diel
integer, allocatable :: wann_bnd(:,:)

if (wannier) then
  allocate(wann_bnd(nstsv,nkptnr))
  wann_bnd=0
! mark all bands that contribute to WF expansion
  do ikloc=1,nkptnr_loc
    ik=iknrglob2(ikloc,mpi_x(1))
    do i=1,nstsv
      do n=1,nwann
        if (abs(wann_c(n,i,ikloc)).gt.1d-10) wann_bnd(i,ik)=1
      enddo
    enddo
  enddo
  call i_reduce_cart(comm_cart_100,.true.,wann_bnd,nstsv*nkptnr)
  do i=1,nstsv
    do ik=1,nkptnr
      wann_bnd(i,1)=wann_bnd(i,1)+wann_bnd(i,ik)
    enddo
  enddo
! find lowest band
  do i=1,nstsv
    if (wann_bnd(i,1).ne.0) then
      min_e1_wann=i 
      exit
    endif
  enddo
! find highest band
  do i=nstsv,1,-1
    if (wann_bnd(i,1).ne.0) then
      max_e2_wann=i 
      exit
    endif
  enddo
  deallocate(wann_bnd)
endif
if (req.and.wproc) then
  write(150,*)
  write(150,'("Bloch functions band interval (N1,N2 or E1,E2) : ",2F8.3)')lr_e1,lr_e2
  if (wannier) then
    write(150,'("Wannier functions band interval (N1,N2 or E1,E2) : ",2F8.3)')min_e1_wann,max_e2_wann
  endif
endif
if (req) then
  nmegqblhmax=0
  min_e12=100.d0
endif
do ikloc=1,nkptnr_loc
  ik=iknrglob2(ikloc,mpi_x(1))
  jk=idxkq(1,ik)
  i=0
  do ist1=1,nstsv
    do ist2=1,nstsv
      le1=bndint(ist1,evalsvnr(ist1,ik),lr_e1,lr_e2)
      le2=bndint(ist2,evalsvnr(ist2,jk),lr_e1,lr_e2)
      d1=occsvnr(ist1,ik)-occsvnr(ist2,jk)
      ldocc=abs(d1).gt.1d-5
      laddme=.false.
! comment:
!   include transition between bands ist1 and ist2 when:
!     1a. we are doing response in Bloch basis and difference of band 
!         occupation numbers in not zero     OR
!     1b. we are doing response in Wannier basis or constrained RPA
!     2.  both bands ist1 and ist2 fall into energy interval
      lwann=.false.
      if (wannier) then
        le1w=bndint(ist1,evalsvnr(ist1,ik),min_e1_wann,max_e2_wann)
        le2w=bndint(ist2,evalsvnr(ist2,jk),min_e1_wann,max_e2_wann)
        if ((lwannresp.and..not.wann_diel()).and.(le1w.and.le2w)) lwann=.true.
        if (crpa.and.(le1w.and.le2w)) lwann=.true.
      endif
      if ((ldocc.or.lwann).and.(le1.and.le2)) then
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
          bmegqblh(1,i,ikloc)=ist1
          bmegqblh(2,i,ikloc)=ist2
          !ime(3,i,ikloc)=1 
          !docc(i,ikloc)=d1
        endif
        if (req) then
          min_e12=min(min_e12,abs(evalsvnr(ist1,ik)-evalsvnr(ist2,jk)))
        endif
      endif
    enddo
  enddo
  if (.not.req) nmegqblh(ikloc)=i
  if (req) nmegqblhmax=max(nmegqblhmax,i)
enddo !ikloc

if (req) then
#ifdef _MPI_
  call d_reduce_cart2(comm_cart_100,.false.,min_e12,1,MPI_MIN)
#endif
  if (wproc) then
    write(150,*)
    write(150,'("Minimal energy transition (eV) : ",F12.6)')min_e12*ha2ev
  endif
endif

return
end

