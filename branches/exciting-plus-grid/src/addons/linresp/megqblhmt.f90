subroutine megqblhmt(nmegqblh_,bmegqblh_,wfsvmt1,wfsvmt2,ngumax,ngu,gu,igu, &
  megqblh_,ik,jk)
use modmain
implicit none
! arguments
integer, intent(in) :: nmegqblh_
integer, intent(in) :: bmegqblh_(2,nmegqblhmax)
integer, intent(in) :: ngumax
integer, intent(in) :: ngu(natmtot,ngvecme)
integer, intent(in) :: igu(4,ngumax,natmtot,ngvecme)
complex(4), intent(in) :: gu(ngumax,natmtot,ngvecme)
complex(8), intent(in) :: wfsvmt1(lmmaxvr,nrfmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfsvmt2(lmmaxvr,nrfmax,natmtot,nspinor,nstsv)
complex(8), intent(inout) :: megqblh_(ngvecme,nmegqblhmax)
integer, intent(in) :: ik
integer, intent(in) :: jk
! local variables
integer ig,i,j,ist1,ist2,ias,io1,io2,lm1,lm2,ispn,ispn2,igloc
integer idx_g1,idx_g2,idx0,bs
complex(8) a1(lmmaxvr,nrfmax,natmtot)
complex(8), allocatable :: megq_tmp(:,:)
logical l1
integer ist1_prev
complex(8), external :: zdotu,zdotc

! TODO: add explanation about how it all works


allocate(megq_tmp(ngvecme,nmegqblhmax))
megq_tmp=zzero

!call idxbos(ngvecme,mpi_dims(2),mpi_x(2)+1,idx0,bs)
!idx_g1=idx0+1
!idx_g2=idx0+bs

ist1_prev=-1
do igloc=1,ngvecmeloc
  ig=mpi_grid_map(ngvecme,dim_g,loc=igloc)
  do ispn=1,nspinor
    ispn2=ispn
    if (lrtype.eq.1) ispn2=3-ispn
    do i=1,nmegqblh_
      ist1=bmegqblh_(1,i)
      ist2=bmegqblh_(2,i)
! precompute
      if (ist1.ne.ist1_prev) then
        a1=zzero
        l1=.true.
        if (spinpol) then
          if (spinor_ud(ispn,ist1,ik).eq.0) l1=.false.
        endif
        if (l1) then
          do ias=1,natmtot
            do j=1,ngu(ias,igloc)
              lm1=igu(1,j,ias,igloc)
              lm2=igu(2,j,ias,igloc)
              io1=igu(3,j,ias,igloc)
              io2=igu(4,j,ias,igloc)
              a1(lm2,io2,ias)=a1(lm2,io2,ias)+&
                dconjg(wfsvmt1(lm1,io1,ias,ispn,ist1))*gu(j,ias,igloc)
            enddo !j
          enddo !ias
        endif
        ist1_prev=ist1
      endif !ist1.ne.ist1_prev
      l1=.true.
      if (spinpol) then
        if (spinor_ud(ispn2,ist2,jk).eq.0) l1=.false.
      endif
      if (l1) then
        do ias=1,natmtot
          megq_tmp(ig,i)=megq_tmp(ig,i)+zdotu(lmmaxvr,a1(1,1,ias),1,&
            wfsvmt2(1,1,ias,ispn2,ist2),1)
          do io2=2,nrfmax
            megq_tmp(ig,i)=megq_tmp(ig,i)+zdotu(16,a1(1,io2,ias),1,&
              wfsvmt2(1,io2,ias,ispn2,ist2),1)    
          enddo
        enddo !ias
      endif
    enddo !i
  enddo !ispn
enddo !igloc

call mpi_grid_reduce(megq_tmp(1,1),ngvecme*nmegqblhmax,dims=(/dim_g/))
megqblh_=megqblh_+megq_tmp

deallocate(megq_tmp)

return
end

