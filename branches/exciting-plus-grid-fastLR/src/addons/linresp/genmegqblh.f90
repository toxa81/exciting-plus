subroutine genmegqblh(ikloc,ngntujumax,ngntuju,igntuju,gntuju,ngknr1,ngknr2,&
  igkignr1,igkignr2,wfsvmt1,wfsvmt2,wfsvit1,wfsvit2)
use modmain
implicit none
integer, intent(in) :: ikloc
integer, intent(in) :: ngntujumax
integer, intent(in) :: ngntuju(nspecies,ngvecme)
integer, intent(in) :: igntuju(4,ngntujumax,nspecies,ngvecme)
integer, intent(in) :: ngknr1
integer, intent(in) :: ngknr2
integer, intent(in) :: igkignr1(ngkmax)
integer, intent(in) :: igkignr2(ngkmax)
complex(8), intent(in) :: gntuju(ngntujumax,nspecies,ngvecme)
complex(8), intent(in) :: wfsvmt1(lmmaxvr,nrfmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfsvmt2(lmmaxvr,nrfmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfsvit1(ngkmax,nspinor,nstsv)
complex(8), intent(in) :: wfsvit2(ngkmax,nspinor,nstsv)

integer wfsize
integer ivg1(3),ivg2(3)
integer i,j,ik,jk,offs,igkq,n1,ispn1,ispn2,ist1,ist2,is
integer ig,ig1,ig2,io1,io2,lm1,lm2,ias
complex(8) zt1
logical l1
logical, allocatable :: lf1g(:)
complex(8), allocatable :: a1(:)
complex(8), allocatable :: wftmp1(:,:)
complex(8), allocatable :: wftmp2(:,:)


wfsize=lmmaxvr*nrfmax*natmtot+ngknr2
allocate(wftmp1(wfsize,ngvecme))
allocate(wftmp2(wfsize,nstsv))
allocate(a1(ngrtot))
allocate(lf1g(ngrtot))

! global k-point
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
! jk=k+q-G_q
jk=idxkq(1,ik)
! G_q vector 
igkq=idxkq(2,ik)
! offset in the global array with interband transitions
offs=nmegqblhloc(2,ikloc)

do ispn1=1,nspinor
  if (lrtype.eq.0) then
    ispn2=ispn1
  else if (lrtype.eq.1) then 
    ispn2=3-ispn1
  endif
! index of local part of interband transitions
  i=1
! go through the local fraction of interband transitions    
  do while (i.le.nmegqblhloc(1,ikloc))
! left <bra| state (referred through local lindex+offset)
    ist1=bmegqblh(1,i+offs,ikloc)
    wftmp1=zzero
    a1=zzero        
    lf1g=.false.
    l1=.true.
    if (spinpol) then
      if (spinor_ud(ispn1,ist1,ik).eq.0) l1=.false.
    endif
    if (l1) then
      do ig=1,ngvecme
! precompute muffint-tin part of \psi_1^{*}(r)*e^{-i(G+q)r}
        do ias=1,natmtot
          is=ias2is(ias)
          do j=1,ngntuju(is,ig)
            lm1=igntuju(1,j,is,ig)
            lm2=igntuju(2,j,is,ig)
            io1=igntuju(3,j,is,ig)
            io2=igntuju(4,j,is,ig)
            wftmp1(lm2+(io2-1)*lmmaxvr+(ias-1)*lmmaxvr*nrfmax,ig)= &
              wftmp1(lm2+(io2-1)*lmmaxvr+(ias-1)*lmmaxvr*nrfmax,ig)+&
              dconjg(wfsvmt1(lm1,io1,ias,ispn1,ist1))*gntuju(j,is,ig)*dconjg(lr_sfacgq0(ig,ias))
          enddo !j
        enddo !ias
! precompute interstitial part of \psi_1^{*}(r)*e^{-i(G+q)r}
        do ig2=1,ngknr2
! ivg2(:)=G+Gq-G2         
          ivg2(:)=ivg(:,ig+gvecme1-1)+ivg(:,igkq)-ivg(:,igkignr2(ig2))
          if (ivg2(1).lt.intgv(1,1).or.ivg2(1).gt.intgv(1,2).or.&
              ivg2(2).lt.intgv(2,1).or.ivg2(2).gt.intgv(2,2).or.&
              ivg2(3).lt.intgv(3,1).or.ivg2(3).gt.intgv(3,2)) then
            write(*,*)
            write(*,'("Error(megqblhit): G-vector is outside of boundaries")')
            write(*,'("  G+G_q-G2 : ",3I5)')ivg2
            write(*,'("  boundaries : ",2I5,",",2I5,",",2I5)')&
              intgv(1,:),intgv(2,:),intgv(3,:)
            write(*,*)
            call pstop
          endif
! if we don't have this Fourier coefficient
          if (.not.lf1g(ivgig(ivg2(1),ivg2(2),ivg2(3)))) then
            zt1=zzero
! compute \sum_{G1} u_1^{*}(G1)*\theta(G1+G)
            do ig1=1,ngknr1
! ivg1(:)=G1+(G+Gq-G2)             
              ivg1(:)=ivg(:,igkignr1(ig1))+ivg2(:)
              if (ivg1(1).lt.intgv(1,1).or.ivg1(1).gt.intgv(1,2).or.&
                  ivg1(2).lt.intgv(2,1).or.ivg1(2).gt.intgv(2,2).or.&
                  ivg1(3).lt.intgv(3,1).or.ivg1(3).gt.intgv(3,2)) then
                write(*,*)
                write(*,'("Error(megqblhit): G-vector is outside of boundaries")')
                write(*,'("  G1+G : ",3I5)')ivg1
                write(*,'("  boundaries : ",2I5,",",2I5,",",2I5)')&
                  intgv(1,:),intgv(2,:),intgv(3,:)
                write(*,*)
                call pstop
              endif
              zt1=zt1+dconjg(wfsvit1(ig1,ispn1,ist1))* &
                      cfunig(ivgig(ivg1(1),ivg1(2),ivg1(3)))
            enddo !ig1
! save the coefficient
            a1(ivgig(ivg2(1),ivg2(2),ivg2(3)))=zt1
! mark it as computed
            lf1g(ivgig(ivg2(1),ivg2(2),ivg2(3)))=.true.
          endif !.not.lf1g(ivgig(ivg2(1),ivg2(2),ivg2(3)))
          wftmp1(lmmaxvr*nrfmax*natmtot+ig2,ig)=a1(ivgig(ivg2(1),ivg2(2),ivg2(3)))
        enddo !ig2
      enddo !ig      
    endif !l1
    n1=0
! collect right |ket> states into matrix wftmp2
    do while ((i+n1).le.nmegqblhloc(1,ikloc))
      if (bmegqblh(1,i+offs+n1,ikloc).ne.bmegqblh(1,i+offs,ikloc)) exit
      ist2=bmegqblh(2,i+offs+n1,ikloc)
      n1=n1+1
      call memcopy(wfsvmt2(1,1,1,ispn2,ist2),wftmp2(1,n1),16*lmmaxvr*nrfmax*natmtot)
      call memcopy(wfsvit2(1,ispn2,ist2),wftmp2(lmmaxvr*nrfmax*natmtot+1,n1),16*ngknr2)
    enddo !while
! update several matrix elements by doing matrix*matrix operation
!  me(ib,ig)=wftmp2(ig2,ib)^{T}*wftmp1(ig2,ig)
    call zgemm('T','N',n1,ngvecme,wfsize,zone,wftmp2,wfsize,wftmp1,wfsize,&
      zone,megqblh(i,1,ikloc),nmegqblhlocmax)
    i=i+n1
  enddo !while
enddo !ispn1
    
deallocate(wftmp1)
deallocate(wftmp2)
deallocate(a1)
deallocate(lf1g)

return
end