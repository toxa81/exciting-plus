subroutine getmeidx(req)
use modmain
use mod_nrkp
implicit none
! arguments
logical, intent(in) :: req

integer i,ik,jk,ist1,ist2,ikloc,n
logical laddme,ldocc
logical l11,l12,l21,l22,le1,le2,lwann,l3,l4
integer, allocatable :: wann_bnd(:,:)
logical, external :: bndint
logical, external :: wann_diel

if (wannier_megq) then
  allocate(wann_bnd(nstsv,nkptnr))
  wann_bnd=0
! mark all bands that contribute to WF expansion
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do i=1,nstsv
      do n=1,nwann
        if (abs(wann_c(n,i,ikloc)).gt.1d-10) wann_bnd(i,ik)=1
      enddo
    enddo
  enddo
  call mpi_grid_reduce(wann_bnd(1,1),nstsv*nkptnr,dims=(/dim_k/),all=.true.,&
    op=op_max)
endif !wannier_megq
if (req) then
  nmegqblhmax=0
  lr_min_e12=100.d0
endif
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  jk=idxkq(1,ik)
  i=0
  do ist1=1,nstsv
    do ist2=1,nstsv
! TODO: more simple code here, currently it is not clear what is going on
!   include transition between bands ist1 and ist2 when:
!     1. difference of band occupation numbers in not zero and
!       both bands ist1 and ist2 fall into energy interval
!     2. this bands are necessary to compute matrix elements
!       in Wannier basis
      le1=bndint(ist1,evalsvnr(ist1,ik),lr_e1,lr_e2)
      le2=bndint(ist2,evalsvnr(ist2,jk),lr_e1,lr_e2)
      ldocc=abs(occsvnr(ist1,ik)-occsvnr(ist2,jk)).gt.1d-6
      lwann=.false.
      if (wannier_megq) then
        l3=(wann_bnd(ist1,ik).ne.0.and.wann_bnd(ist2,jk).ne.0)
        if (.not.wann_diel().and.l3) lwann=.true.
        if (wann_diel().and.ldocc.and.l3) lwann=.true.
        if (all_wan_ibt.and.l3) lwann=.true.
      endif
      l4=.false.
! for response and cRPA we need bands in [e1,e2] interval 
      if (task.eq.400.or.task.eq.401) then
        if (ldocc.and.le1.and.le2) l4=.true.
      endif
! bands for matrix in Wannier basis
      if (lwann) l4=.true.
      laddme=.false.
      if (l4) then
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
          if (wannier_megq) then
            if (wann_bnd(ist1,ik).ne.0.and.wann_bnd(ist2,jk).ne.0) then
              nmegqblhwan(ikloc)=nmegqblhwan(ikloc)+1
              imegqblhwan(nmegqblhwan(ikloc),ikloc)=i
            endif
          endif
        endif
        if (req) then
          lr_min_e12=min(lr_min_e12,abs(evalsvnr(ist1,ik)-evalsvnr(ist2,jk)))
        endif
      endif
    enddo
  enddo
  if (.not.req) nmegqblh(ikloc)=i
  if (req) nmegqblhmax=max(nmegqblhmax,i)
enddo !ikloc

if (req) then
  call mpi_grid_reduce(lr_min_e12,dims=(/dim_k/),op=op_min)
endif

if (wannier_megq) then
  deallocate(wann_bnd)
endif

return
end

