subroutine genmegqblh(iq,ikloc,ngknr1,ngknr2,igkignr1,igkignr2,wfsvmt1,wfsvmt2,&
  wfsvit1,wfsvit2)
use modmain
use mod_addons_q
implicit none
integer, intent(in) :: iq
integer, intent(in) :: ikloc
integer, intent(in) :: ngknr1
integer, intent(in) :: ngknr2
integer, intent(in) :: igkignr1(ngkmax)
integer, intent(in) :: igkignr2(ngkmax)
complex(8), intent(in) :: wfsvmt1(lmmaxvr,nufrmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfsvmt2(lmmaxvr,nufrmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfsvit1(ngkmax,nspinor,nstsv)
complex(8), intent(in) :: wfsvit2(ngkmax,nspinor,nstsv)

integer wfsize
integer ivg1(3)
integer i,j,ik,jk,offs,igkq,n1,ispn1,ispn2,ist1,ist2,ic
integer ig,ig1,ig2,io1,io2,lm1,lm2,ias,ifg,ir
logical l1
logical, allocatable :: lf1g(:)
complex(8), allocatable :: a1(:)
complex(8), allocatable :: wftmp1(:,:)
complex(8), allocatable :: wftmp2(:,:)
complex(8) b1(lmmaxvr,nufrmax)
complex(8), allocatable :: wfir1(:)


wfsize=lmmaxvr*nufrmax*natmtot+ngknr2
allocate(wftmp1(wfsize,ngvecme))
allocate(wftmp2(wfsize,nstsv))
allocate(a1(ngrtot))
allocate(lf1g(ngrtot))
allocate(wfir1(ngrtot))

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
        call timer_start(3)
        do ias=1,natmtot
          b1=dconjg(wfsvmt1(:,:,ias,ispn1,ist1)*sfacgq(ig,ias))
          ic=ias2ic(ias)
          do j=1,ngntuju(ic,ig)
            lm1=igntuju(1,j,ic,ig)
            lm2=igntuju(2,j,ic,ig)
            io1=igntuju(3,j,ic,ig)
            io2=igntuju(4,j,ic,ig)
            wftmp1(lm2+(io2-1)*lmmaxvr+(ias-1)*lmmaxvr*nufrmax,ig)= &
              wftmp1(lm2+(io2-1)*lmmaxvr+(ias-1)*lmmaxvr*nufrmax,ig)+&
              b1(lm1,io1)*gntuju(j,ic,ig)
          enddo !j
        enddo !ias
        call timer_stop(3)
      enddo !ig  
! interstitial part
      call timer_start(4)
      wfir1=zzero
      do ig1=1,ngknr1
        ifg=igfft(igkignr1(ig1))
        wfir1(ifg)=dconjg(wfsvit1(ig1,ispn1,ist1))
      enddo
      call zfftifc(3,ngrid,1,wfir1)
      do ir=1,ngrtot
        wfir1(ir)=wfir1(ir)*cfunir(ir)
      enddo
      call zfftifc(3,ngrid,-1,wfir1)
      do ig=1,ngvecme
        do ig2=1,ngknr2
          ivg1(:)=ivg(:,igkignr2(ig2))-ivg(:,igqig(ig,iq))-ivg(:,igkq)
          ifg=igfft(ivgig(ivg1(1),ivg1(2),ivg1(3)))
          wftmp1(lmmaxvr*nufrmax*natmtot+ig2,ig)=wfir1(ifg)
        enddo
      enddo
      call timer_stop(4)      
    endif !l1
    call timer_start(5)
    n1=0
! collect right |ket> states into matrix wftmp2
    do while ((i+n1).le.nmegqblhloc(1,ikloc))
      if (bmegqblh(1,i+offs+n1,ikloc).ne.bmegqblh(1,i+offs,ikloc)) exit
      ist2=bmegqblh(2,i+offs+n1,ikloc)
      n1=n1+1
      call memcopy(wfsvmt2(1,1,1,ispn2,ist2),wftmp2(1,n1),16*lmmaxvr*nufrmax*natmtot)
      call memcopy(wfsvit2(1,ispn2,ist2),wftmp2(lmmaxvr*nufrmax*natmtot+1,n1),16*ngknr2)
    enddo !while
! update several matrix elements by doing matrix*matrix operation
!  me(ib,ig)=wftmp2(ig2,ib)^{T}*wftmp1(ig2,ig)
    call zgemm('T','N',n1,ngvecme,wfsize,zone,wftmp2,wfsize,wftmp1,wfsize,&
      zone,megqblh(i,1,ikloc),nmegqblhlocmax)
    i=i+n1
    call timer_stop(5)
  enddo !while
enddo !ispn1
    
deallocate(wftmp1)
deallocate(wftmp2)
deallocate(a1)
deallocate(lf1g)
deallocate(wfir1)

return
end