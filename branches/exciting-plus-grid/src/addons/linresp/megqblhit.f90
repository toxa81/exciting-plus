subroutine megqblhit(nmegqblh_,bmegqblh_,ngknr1,ngknr2,igkignr1,igkignr2, &
  igkq,wfsvit1,wfsvit2,ik,jk,megqblh_)
use modmain
implicit none
! arguments
integer, intent(in) :: nmegqblh_
integer, intent(in) :: bmegqblh_(2,nmegqblhmax)
integer, intent(in) :: ngknr1
integer, intent(in) :: ngknr2
integer, intent(in) :: igkignr1(ngkmax)
integer, intent(in) :: igkignr2(ngkmax)
integer, intent(in) :: igkq
complex(8), intent(in) :: wfsvit1(ngkmax,nspinor,nstsv)
complex(8), intent(in) :: wfsvit2(ngkmax,nspinor,nstsv)
integer, intent(in) :: ik
integer, intent(in) :: jk
complex(8), intent(inout) :: megqblh_(ngvecme,nmegqblhmax)

complex(8), allocatable :: mit(:,:)
complex(8), allocatable :: a1(:) 
complex(8), allocatable :: megq_tmp(:,:)
integer ivg1(3),ig1,ig2
integer ig,ist1,ist2,i,ispn,ispn2,j,ist
integer idx0,bs,idx_g1,idx_g2,igp,ifg,ir,idx_i1,idx_i2
complex(8) zt1
logical l1
complex(8), external :: zdotu,zdotc
integer ist1_prev
logical, allocatable :: lig1(:)

! TODO: add explanation about how it all works

allocate(megq_tmp(ngvecme,nmegqblhmax))
megq_tmp=zzero

allocate(a1(ngrtot))
allocate(lig1(ngrtot))
! split interband transitions between processors
bs=mpi_grid_map(nmegqblh_,dim2,offs=idx0)
idx_i1=idx0+1
idx_i2=idx0+bs
! no previos state was computed  
ist1_prev=-1
do ispn=1,nspinor
  ispn2=ispn
  if (lrtype.eq.1) ispn2=3-ispn
  do i=idx_i1,idx_i2
    ist1=bmegqblh_(1,i)
    ist2=bmegqblh_(2,i)
! precompute
    if (ist1.ne.ist1_prev) then
      a1=zzero
      lig1=.false.
      l1=.true.
      if (spinpol) then
        if (spinor_ud(ispn,ist1,ik).eq.0) l1=.false.
      endif
      if (l1) then
! determine, which Fourie coefficients we need for the product u_1^{*}(r)*\theta(r)
        do ig=1,ngvecme
          do ig2=1,ngknr2
! G1=G+Gq-G2
            ivg1(:)=ivg(:,ig+gvecme1-1)+ivg(:,igkq)-ivg(:,igkignr2(ig2))
            if (ivg1(1).lt.intgv(1,1).or.ivg1(1).gt.intgv(1,2).or.&
                ivg1(2).lt.intgv(2,1).or.ivg1(1).gt.intgv(2,2).or.&
                ivg1(3).lt.intgv(3,1).or.ivg1(1).gt.intgv(3,2)) then
              write(*,*)
              write(*,'("Error(megqblhit): G-vector is outside of boundaries")')
              write(*,'("  G+G_q-G2 : ",3I5)')ivg1
              write(*,'("  boundaries : ",2I5,",",2I5,",",2I5)')&
                intgv(1,:),intgv(2,:),intgv(3,:)
              write(*,*)
              call pstop
            endif
            lig1(ivgig(ivg1(1),ivg1(2),ivg1(3)))=.true.
          enddo
        enddo !ig 
        do ig=1,ngrtot
          if (lig1(ig)) then
            zt1=zzero
            do ig1=1,ngknr1
! Gt=G1+G
              ivg1(:)=ivg(:,ig)+ivg(:,igkignr1(ig1))
              zt1=zt1+dconjg(wfsvit1(ig1,ispn,ist1))* &
                      cfunig(ivgig(ivg1(1),ivg1(2),ivg1(3)))
            enddo
            a1(ig)=zt1
          endif
        enddo
      endif
      ist1_prev=ist1
    endif !ist1.ne.ist1_prev
    l1=.true.
    if (spinpol) then
      if (spinor_ud(ispn2,ist2,jk).eq.0) l1=.false.
    endif
    if (l1) then
      do ig=1,ngvecme
        zt1=zzero
        do ig2=1,ngknr2
          ivg1(:)=ivg(:,ig+gvecme1-1)+ivg(:,igkq)-ivg(:,igkignr2(ig2))
          zt1=zt1+a1(ivgig(ivg1(1),ivg1(2),ivg1(3)))*wfsvit2(ig2,ispn2,ist2)
        enddo !ig2
        megq_tmp(ig,i)=megq_tmp(ig,i)+zt1
      enddo !ig
    endif
  enddo !i
enddo !ispn
deallocate(a1,lig1)
if (mpi_grid_size(dim2).gt.1.) then
  call mpi_grid_reduce(megq_tmp(1,1),ngvecme*nmegqblhmax,dims=(/dim2/))
endif
megqblh_=megqblh_+megq_tmp
deallocate(megq_tmp)

return
end