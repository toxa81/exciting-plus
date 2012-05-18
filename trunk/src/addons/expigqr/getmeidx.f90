subroutine getmeidx
use modmain
use mod_nrkp
use mod_expigqr
implicit none
! local variables
integer i,ik,jk,ist1,ist2,ikloc,n,i1,i2,n1,n2
logical laddibt,ldocc,laddme,lwanibt
logical l11,l12,l21,l22,le1,le2
integer*2, allocatable :: wann_bnd_n(:,:)
integer*2, allocatable :: wann_bnd_k(:,:)
logical, external :: bndint

if (wannier_megq) then
  allocate(wann_bnd_n(nstsv,nwantot))
  allocate(wann_bnd_k(nstsv,nkptnr))
  wann_bnd_n=0
  wann_bnd_k=0
! mark all bands that contribute to WF expansion
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do i=1,nstsv
      do n=1,nwantot
        if (abs(wanncnrloc(n,i,ikloc)).gt.1d-10) then
! at least at one k-point the state i contributes to WF n
          wann_bnd_n(i,n)=1
! for k-point ik band i contributes at least to one WF 
          wann_bnd_k(i,ik)=1
        endif
      enddo
    enddo
  enddo
  call mpi_grid_reduce(wann_bnd_n(1,1),nstsv*nwantot,dims=(/dim_k/),&
    &all=.true.,op=op_max)
  call mpi_grid_reduce(wann_bnd_k(1,1),nstsv*nkptnr,dims=(/dim_k/),&
    &all=.true.,op=op_max)
endif !wannier_megq
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  jk=idxkq(1,ik)
  i=0
  do ist1=1,nstsv
    do ist2=1,nstsv
      lwanibt=.false.
! include transition between bands ist1 and ist2 when:
!  1) both bands ist1 and ist2 fall into energy interval and
!     difference of band occupation numbers in not zero  
      le1=bndint(ist1,evalsvnr(ist1,ik),megq_include_bands(1),&
        megq_include_bands(2))
      le2=bndint(ist2,evalsvnr(ist2,jk),megq_include_bands(1),&
        megq_include_bands(2))
      ldocc=abs(occsvnr(ist1,ik)-occsvnr(ist2,jk)).gt.1d-6
      laddibt=(le1.and.le2.and.ldocc)
!  2) this bands are necessary to compute matrix elements in Wannier basis
      if (wannier_megq) then
! if this bands contribute to at least one Wannier function
        if (wann_bnd_k(ist1,ik).ne.0.and.wann_bnd_k(ist2,jk).ne.0) then
! check contribution to each Wannier function
! NOTE: current implementation works with small overhead; in principle we
!  must check Bloch contribution to particular WF at particular k-point, but
!  this requires much bigger array of "Bloch contribution to WF" flags
          do i1=1,megqwantran%nwan
            n1=megqwantran%iwan(i1)
            do i2=1,megqwantran%nwan
              n2=megqwantran%iwan(i2)
              if (megqwantran%wt(n1,n2).eq.1.and.wann_bnd_n(ist1,n1).eq.1.and.&
                  &wann_bnd_n(ist2,n2).eq.1) then
                laddibt=.true.
                lwanibt=.true.
              endif
            enddo !i2
          enddo !i1
        endif
      endif !wannier_megq
      laddme=.false.
! final check: don't add matrix element if it is zero
      if (laddibt) then
        if (.not.spinpol) then
          laddme=.true.
        else
          l11=spinor_ud(1,ist1,ik).eq.1.and.spinor_ud(1,ist2,jk).eq.1
          l12=spinor_ud(1,ist1,ik).eq.1.and.spinor_ud(2,ist2,jk).eq.1
          l21=spinor_ud(2,ist1,ik).eq.1.and.spinor_ud(1,ist2,jk).eq.1
          l22=spinor_ud(2,ist1,ik).eq.1.and.spinor_ud(2,ist2,jk).eq.1
          if (expigqr22.eq.1.and.(l11.or.l22)) laddme=.true.
          !if (expigqr22.eq.2.and.(l12.or.l21)) laddme=.true.
        endif
      endif
      if (laddme) then
        i=i+1
        bmegqblh(1,i,ikloc)=ist1
        bmegqblh(2,i,ikloc)=ist2
        if (lwanibt) then
          nmegqblhwan(ikloc)=nmegqblhwan(ikloc)+1
          imegqblhwan(nmegqblhwan(ikloc),ikloc)=i
        endif
      endif !laddme
    enddo !ist2
  enddo !ist1
  nmegqblh(ikloc)=i
enddo !ikloc
if (wannier_megq) then
  deallocate(wann_bnd_n)
  deallocate(wann_bnd_k)
endif
return
end

