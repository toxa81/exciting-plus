subroutine megqblhit(nmegqblh_,bmegqblh_,ngknr1,ngknr2,igkignr1,igkignr2, &
  igkq,wfsvit1,wfsvit2,megqblh_,gvit,ik,jk)
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
complex(8), intent(inout) :: megqblh_(ngvecme,nmegqblhmax)
complex(8), intent(in) :: gvit(intgv(1,1):intgv(1,2),intgv(2,1):intgv(2,2),&
                            intgv(3,1):intgv(3,2))
integer, intent(in) :: ik
integer, intent(in) :: jk

complex(8), allocatable :: mit(:,:)
complex(8), allocatable :: a1(:) 
complex(8), allocatable :: megq_tmp(:,:)
complex(8), allocatable :: wfir1(:)
integer ivg1(3),ig1,ig2
integer ig,ist1,ist2,i,ispn,ispn2,j,ist
integer idx0,bs,idx_g1,idx_g2,igp,ifg,ir,idx_i1,idx_i2
complex(8) zt1
logical l1
complex(8), external :: zdotu,zdotc
integer ist1_prev

! TODO: add explanation about how it all works

allocate(megq_tmp(ngvecme,nmegqblhmax))
megq_tmp=zzero

! analytical expression for interstitial contributuion : fast for small
!  number of G-vectors for wave-functions and matrix elements (typically
!  this is the case of small systems like Al or NiO)
if (.not.lfftit) then
  allocate(a1(ngknr2))
  allocate(mit(ngknr1,ngknr2))
! split G-vectors between processors
  call idxbos(ngvecme,mpi_dims(2),mpi_x(2)+1,idx0,bs)
  idx_g1=idx0+1
  idx_g2=idx0+bs
! no previos state was computed  
  ist1_prev=-1
  do ig=idx_g1,idx_g2
	call genpwit2(ngknr1,ngknr2,igkignr1,igkignr2,-(ivg(:,ig+gvecme1-1)+&
	  ivg(:,igkq)),gvit,mit)
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
			call zgemv('T',ngknr1,ngknr2,zone,mit,ngknr1,&
			  dconjg(wfsvit1(:,ispn,ist1)),1,zzero,a1,1)
		  endif
		  ist1_prev=ist1
		endif
		l1=.true.
		if (spinpol) then
		  if (spinor_ud(ispn2,ist2,jk).eq.0) l1=.false.
		endif
		if (l1) then
		  megq_tmp(ig,i)=megq_tmp(ig,i)+zdotu(ngknr2,wfsvit2(1,ispn2,ist2),1,a1,1)
		endif
	  enddo !i  
	enddo !ispn
  enddo !ig
  deallocate(mit,a1)
! numerical expression for interstitial contributuion using the FFT technique
!   this is much faster for the case of "real" systems with lots of atoms
!   and G-vectors for matrix elements; slow for small systems due to FFT
!   overhead
else
  allocate(wfir1(ngrtot))
! split interband transitions between processors
  call idxbos(nmegqblh_,mpi_dims(2),mpi_x(2)+1,idx0,bs)
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
		wfir1=zzero
		l1=.true.
		if (spinpol) then
		  if (spinor_ud(ispn,ist1,ik).eq.0) l1=.false.
		endif
		if (l1) then
! compute u_{nk}(r)
          do ig1=1,ngknr1
			ifg=igfft(igkignr1(ig1))
			wfir1(ifg)=wfsvit1(ig1,ispn,ist1)
		  enddo
		  call zfftifc(3,ngrid,1,wfir1)
! compute f_{nk}(r)=u_{nk}(r)*\theta(r)
		  do ir=1,ngrtot
			wfir1(ir)=wfir1(ir)*cfunir(ir)
		  enddo
! compute f_{nk}(G)
		  call zfftifc(3,ngrid,-1,wfir1)
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
			ivg1(:)=ivg(:,igkignr2(ig2))-ivg(:,ig)-ivg(:,igkq)
			ifg=igfft(ivgig(ivg1(1),ivg1(2),ivg1(3)))
			zt1=zt1+dconjg(wfir1(ifg))*wfsvit2(ig2,ispn2,ist2)
		  enddo !ig2
		  megq_tmp(ig,i)=megq_tmp(ig,i)+zt1
		enddo !ig
	  endif
	enddo !i
  enddo !ispn
  deallocate(wfir1)
endif

if (mpi_dims(2).gt.1) then
  call d_reduce_cart(comm_cart_010,.false.,megq_tmp,2*ngvecme*nmegqblhmax)
endif
megqblh_=megqblh_+megq_tmp
deallocate(megq_tmp)

return
end