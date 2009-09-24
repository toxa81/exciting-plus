subroutine zrhoftmt(nme0,ime0,wfsvmt1,wfsvmt2,ngumax,ngu,gu,igu, &
  zrhofc0,ik,jk)
use modmain
implicit none
! arguments
integer, intent(in) :: nme0
integer, intent(in) :: ime0(3,nmemax)
integer, intent(in) :: ngumax
integer, intent(in) :: ngu(natmtot,ngvecme)
integer, intent(in) :: igu(4,ngumax,natmtot,ngvecme)
complex(4), intent(in) :: gu(ngumax,natmtot,ngvecme)
complex(8), intent(in) :: wfsvmt1(lmmaxvr,nrfmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfsvmt2(lmmaxvr,nrfmax,natmtot,nspinor,nstsv)
complex(8), intent(inout) :: zrhofc0(ngvecme,nmemax)
integer, intent(in) :: ik
integer, intent(in) :: jk
! local variables
integer ig,i,j,ist1,ist2,ias,io1,io2,lm1,lm2,ispn,ispn2
integer idx_g1,idx_g2,idx0,bs
complex(8) a1(lmmaxvr,nrfmax),a2(lmmaxvr,nrfmax)
complex(8), allocatable :: zrhofc_tmp(:,:)
logical l1


allocate(zrhofc_tmp(ngvecme,nmemax))
zrhofc_tmp=dcmplx(0.d0,0.d0)

call idxbos(ngvecme,mpi_dims(2),mpi_x(2)+1,idx0,bs)
idx_g1=idx0+1
idx_g2=idx0+bs

do ig=idx_g1,idx_g2
  do ispn=1,nspinor
    ispn2=ispn
    if (lrtype.eq.1) then
      ispn2=3-ispn
    endif
    do i=1,nme0
      ist1=ime0(1,i)
      ist2=ime0(2,i)
      l1=.true.
      if (spinpol) then
        if (spinor_ud(ispn,ist1,ik).eq.0.or.spinor_ud(ispn2,ist2,jk).eq.0) l1=.false.
      endif
      if (l1) then
        do ias=1,natmtot
          a1=dconjg(wfsvmt1(:,:,ias,ispn,ist1))
          a2=wfsvmt2(:,:,ias,ispn2,ist2)
          do j=1,ngu(ias,ig)
            lm1=igu(1,j,ias,ig)
            lm2=igu(2,j,ias,ig)
            io1=igu(3,j,ias,ig)
            io2=igu(4,j,ias,ig)
            zrhofc_tmp(ig,i)=zrhofc_tmp(ig,i)+a1(lm1,io1)*a2(lm2,io2)*gu(j,ias,ig)
          enddo
        enddo !ias
      endif
    enddo !i
  enddo
enddo !ig    

if (mpi_dims(2).gt.1) then
  call d_reduce_cart(comm_cart_010,.false.,zrhofc_tmp,2*ngvecme*nmemax)
endif
zrhofc0=zrhofc0+zrhofc_tmp

deallocate(zrhofc_tmp)

return
end

