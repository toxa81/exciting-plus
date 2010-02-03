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
integer idx_g1,idx_g2,idx0,bs
complex(8) a1(lmmaxvr,nrfmax,natmtot)
complex(8), allocatable :: megq_tmp(:,:)
complex(8), allocatable :: megq_tmp2(:,:)
complex(8), allocatable :: wftmp1(:,:,:)
complex(8), allocatable :: wftmp2(:,:,:,:)
logical l1
integer ist1_prev,n1,n2,ik,jk
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

ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
jk=idxkq(1,ik)


allocate(megq_tmp(ngvecme,nmegqblhmax))
allocate(megq_tmp2(nmegqblhmax,ngvecme))
allocate(wftmp1(lmmaxvr,nrfmax,natmtot))
allocate(wftmp2(lmmaxvr,nrfmax,natmtot,nstsv))

megq_tmp=zzero
megq_tmp2=zzero

ist1_prev=-1
do ig=1,ngvecme
  do ispn1=1,nspinor
    if (lrtype.eq.0) then
      ispn2=ispn1
    else if (lrtype.eq.1) then 
      ispn2=3-ispn1
    endif
    
    
    
    i=1 
    do while (i.le.nmegqblh_)
! precompute
      ist1=bmegqblh_(1,i)
      wftmp1=zzero
      l1=.true.
      if (spinpol) then
        if (spinor_ud(ispn1,ist1,ik).eq.0) l1=.false.
      endif
      if (l1) then
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
      endif
      n1=0
      do while (bmegqblh_(1,i+n1).eq.bmegqblh_(1,i))
        ist2=bmegqblh_(2,i+n1)
        n1=n1+1
        wftmp2(:,:,:,n1)=wfsvmt2(:,:,:,ispn2,ist2)
      enddo
      call zgemv('T',lmmaxvr*nrfmax*natmtot,n1,zone,wftmp2,lmmaxvr*nrfmax*natmtot,&
        wftmp1,1,zone,megq_tmp2(i,ig),1)
      i=i+n1
    enddo
!    
!    
!    
!    do i=1,nmegqblh_
!      ist1=bmegqblh_(1,i)
!      ist2=bmegqblh_(2,i)
!! precompute
!      if (ist1.ne.ist1_prev) then
!        a1=zzero
!        l1=.true.
!        if (spinpol) then
!          if (spinor_ud(ispn1,ist1,ik).eq.0) l1=.false.
!        endif
!        if (l1) then
!          do ias=1,natmtot
!            do j=1,ngntuju(ias,igloc)
!              lm1=igntuju(1,j,ias,igloc)
!              lm2=igntuju(2,j,ias,igloc)
!              io1=igntuju(3,j,ias,igloc)
!              io2=igntuju(4,j,ias,igloc)
!              a1(lm2,io2,ias)=a1(lm2,io2,ias)+&
!                dconjg(wfsvmt1(lm1,io1,ias,ispn1,ist1))*gntuju(j,ias,igloc)
!            enddo !j
!          enddo !ias
!        endif
!        ist1_prev=ist1
!
!      endif !ist1.ne.ist1_prev
!      l1=.true.
!      if (spinpol) then
!        if (spinor_ud(ispn2,ist2,jk).eq.0) l1=.false.
!      endif
!      if (l1) then
!        do ias=1,natmtot
!          do io2=1,nrfmax
!          megq_tmp(ig,i)=megq_tmp(ig,i)+zdotu(lmmaxvr,a1(1,io2,ias),1,&
!            wfsvmt2(1,io2,ias,ispn2,ist2),1)
!          enddo
!          !do io2=2,nrfmax
!          !  megq_tmp(ig,i)=megq_tmp(ig,i)+zdotu(16,a1(1,io2,ias),1,&
!          !    wfsvmt2(1,io2,ias,ispn2,ist2),1)    
!          !enddo
!        enddo !ias
!      endif
!    enddo !i
  enddo !ispn1
  
  !write(*,*)megq_tmp(ig,:)-megq_tmp2(:,ig) 
  megq_tmp(ig,:)=megq_tmp2(:,ig)
enddo !ig

call mpi_grid_reduce(megq_tmp(1,1),ngvecme*nmegqblhmax,dims=(/dim_g/))
megqblh_=megqblh_+megq_tmp

deallocate(megq_tmp,megq_tmp2,wftmp1,wftmp2)

return
end

