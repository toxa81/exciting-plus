subroutine megqblhmt(ikloc,wfsvmt1,wfsvmt2,ngntujumax,ngntuju,igntuju,gntuju)
use modmain
implicit none
! arguments
integer, intent(in) :: ikloc
complex(8), intent(in) :: wfsvmt1(lmmaxvr,nrfmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfsvmt2(lmmaxvr,nrfmax,natmtot,nspinor,nstsv)
integer, intent(in) :: ngntujumax
integer, intent(in) :: ngntuju(natmtot,ngvecme)
integer, intent(in) :: igntuju(4,ngntujumax,natmtot,ngvecme)
complex(8), intent(in) :: gntuju(ngntujumax,natmtot,ngvecme)
! local variables
integer ig,i,j,ist1,ist2,ias,io1,io2,lm1,lm2,ispn1,ispn2
complex(8), allocatable :: wftmp1(:,:,:)
complex(8), allocatable :: wftmp2(:,:,:,:)
logical l1,l2
integer ist1_prev,n1,n2,ik,jk,offs
complex(8), external :: zdotu,zdotc

! TODO: add explanation about how it all works
!write(*,*)'ik=',ik
!write(*,*)'nmegqblh_=',nmegqblh_
!write(*,*)'bmegqblh_=',sum(bmegqblh_)
!write(*,*)'wfsvmt1=',sum(wfsvmt1)
!write(*,*)'wfsvmt2=',sum(wfsvmt2)
!write(*,*)'ngntujumax=',ngntujumax
!write(*,*)'ngntuju=',sum(ngntuju)
!write(*,*)'igntuju=',sum(igntuju)
!write(*,*)'gntuju=',sum(gntuju)
!write(*,*)

! global k-point
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
! jk=k+q-G_q
jk=idxkq(1,ik)
! offset in the global array with interband transitions
offs=nmegqblhloc(2,ikloc)
allocate(wftmp1(lmmaxvr,nrfmax,natmtot))
allocate(wftmp2(lmmaxvr,nrfmax,natmtot,nstsv))
do ig=1,ngvecme
  do ispn1=1,nspinor
    if (lrtype.eq.0) then
      ispn2=ispn1
    else if (lrtype.eq.1) then 
      ispn2=3-ispn1
    endif
! index of local part of interband transitions
    i=1
! go through the own fraction of interband transitions    
    do while (i.le.nmegqblhloc(1,ikloc))
! left <bra| state (referred through local lindex+offset)
      ist1=bmegqblh(1,i+offs,ikloc)
      wftmp1=zzero
      l1=.true.
      if (spinpol) then
        if (spinor_ud(ispn1,ist1,ik).eq.0) l1=.false.
      endif
      if (l1) then
! precompute \psi_1^{*}(r)*e^{-i(G+q)r}
        do ias=1,natmtot
          do j=1,ngntuju(ias,ig)
            lm1=igntuju(1,j,ias,ig)
            lm2=igntuju(2,j,ias,ig)
            io1=igntuju(3,j,ias,ig)
            io2=igntuju(4,j,ias,ig)
            wftmp1(lm2,io2,ias)=wftmp1(lm2,io2,ias)+&
              dconjg(wfsvmt1(lm1,io1,ias,ispn1,ist1))*gntuju(j,ias,ig)
          enddo !j
        enddo !ias
      endif !l1
      n1=0
      l2=.true.
! collect right |ket> states into matrix wftmp2
      do while (l2.and.(i+n1).le.nmegqblhloc(1,ikloc))
        if (bmegqblh(1,i+offs+n1,ikloc).ne.bmegqblh(1,i+offs,ikloc)) l2=.false.
        ist2=bmegqblh(2,i+offs+n1,ikloc)
        n1=n1+1
        wftmp2(:,:,:,n1)=wfsvmt2(:,:,:,ispn2,ist2)
      enddo
! update several matrix elements by doing matrix*vector operation
      call zgemv('T',lmmaxvr*nrfmax*natmtot,n1,zone,wftmp2,lmmaxvr*nrfmax*natmtot,&
        wftmp1,1,zone,megqblh(i,ig,ikloc),1)
      i=i+n1
    enddo
  enddo !ispn1
enddo !ig
deallocate(wftmp1,wftmp2)
return
end

