subroutine zrhoftit(nme0,ime0,ngknr1,ngknr2,igkignr1,igkignr2, &
  igkq,wfsvit1,wfsvit2,zrhofc0,evecfv1,evecsv1,evecfv2,evecsv2,igfft1)
use modmain
implicit none
! arguments
integer, intent(in) :: nme0
integer, intent(in) :: ime0(3,nmemax)
integer, intent(in) :: ngknr1
integer, intent(in) :: ngknr2
integer, intent(in) :: igkq
integer, intent(in) :: igkignr1(ngkmax)
integer, intent(in) :: igkignr2(ngkmax)
complex(8), intent(in) :: wfsvit1(ngkmax,nstsv,nspinor)
complex(8), intent(in) :: wfsvit2(ngkmax,nstsv,nspinor)
complex(8), intent(inout) :: zrhofc0(ngvecme,nmemax)
complex(8), intent(in) :: evecfv1(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: evecsv1(nstsv,nstsv)
complex(8), intent(in) :: evecfv2(nmatmax,nstfv,nspnfv)
complex(8), intent(in) :: evecsv2(nstsv,nstsv)
integer, intent(in) :: igfft1(ngvecme)


complex(8), allocatable :: mit(:,:)
complex(8), allocatable :: a(:,:,:) 
complex(8), allocatable :: zrhofc_tmp(:,:)
complex(8), allocatable :: wfir1(:,:,:)
complex(8), allocatable :: wfir2(:,:,:)
complex(8), allocatable :: zrhoir(:)


integer is,ia,ias,ig,ig1,ig2,ist1,ist2,i,ispn,ispn2,j,ist
integer idx0,bs,idx_g1,idx_g2,igp,ifg,ir,i1,i2
integer iv3g(3)
real(8) v1(3),v2(3),tp3g(2),len3g,t1
complex(8) sfac3g(natmtot),zt1



if (.false.) then
allocate(a(ngknr2,nstsv,nspinor))
allocate(zrhofc_tmp(ngvecme,nmemax))
zrhofc_tmp=dcmplx(0.d0,0.d0)

call idxbos(ngvecme,mpi_dims(2),mpi_x(2)+1,idx0,bs)
idx_g1=idx0+1
idx_g2=idx0+bs

allocate(mit(ngknr1,ngknr2))

do ig=idx_g1,idx_g2
  mit=dcmplx(0.d0,0.d0)
  do ig1=1,ngknr1
    do ig2=1,ngknr2
      ! G1-G2+G+K
      iv3g(:)=ivg(:,igkignr1(ig1))-ivg(:,igkignr2(ig2))+ivg(:,ig+gvecme1-1)+ivg(:,igkq)
      if (sum(abs(iv3g)).eq.0) mit(ig1,ig2)=dcmplx(1.d0,0.d0)
      v2(:)=1.d0*iv3g(:)
      call r3mv(bvec,v2,v1)
      call sphcrd(v1,len3g,tp3g)
      call gensfacgp(1,v1,1,sfac3g)
      do is=1,nspecies
        do ia=1,natoms(is)
	  ias=idxas(ia,is)
	  if (len3g.lt.1d-8) then
	    mit(ig1,ig2)=mit(ig1,ig2)-(fourpi/omega)*dconjg(sfac3g(ias))*(rmt(is)**3)/3.d0
	  else
	    mit(ig1,ig2)=mit(ig1,ig2)-(fourpi/omega)*dconjg(sfac3g(ias)) * &
	      (-(rmt(is)/len3g**2)*cos(len3g*rmt(is))+(1/len3g**3)*sin(len3g*rmt(is)))
	  endif
	enddo !ia
      enddo !is
    enddo
  enddo
  
  a=dcmplx(0.d0,0.d0)
  do ispn=1,nspinor
    do i=1,nstsv
      do ig2=1,ngknr2
        do ig1=1,ngknr1
          a(ig2,i,ispn)=a(ig2,i,ispn) + &
	    dconjg(wfsvit1(ig1,i,ispn))*mit(ig1,ig2)
        enddo
      enddo
    enddo
  enddo
  
  do ispn=1,nspinor
    if (lrtype.eq.0) then
      ispn2=ispn
    endif
    if (lrtype.eq.1) then
      ispn2=3-ispn
    endif
    do i=1,nme0
      ist1=ime0(1,i)
      ist2=ime0(2,i)
      do ig2=1,ngknr2
        zrhofc_tmp(ig,i)=zrhofc_tmp(ig,i)+wfsvit2(ig2,ist2,ispn2)*a(ig2,ist1,ispn)
      enddo
    enddo
  enddo
enddo !ig

if (mpi_dims(2).gt.1) then
  call d_reduce_cart(comm_cart_010,.false.,zrhofc_tmp,2*ngvecme*nmemax)
endif
zrhofc0=zrhofc0+zrhofc_tmp

deallocate(mit,a,zrhofc_tmp)



else

allocate(wfir1(ngrtot,nspinor,nstsv))
allocate(wfir2(ngrtot,nspinor,nstsv))
allocate(zrhoir(ngrtot))
!-----------------------------------!
!     interstitial wavefunction     !
!-----------------------------------!
!t1=1.d0 /sqrt(omega)
wfir1(:,:,:)=dcmplx(0.d0,0.d0)
wfir2(:,:,:)=dcmplx(0.d0,0.d0)

do j=1,nstsv
  if (tevecsv) then
! generate spinor wavefunction from second-variational eigenvectors
    i=0
    do ispn=1,nspinor
      do ist=1,nstfv
        i=i+1
        zt1=evecsv1(i,j)
        if (abs(dble(zt1))+abs(aimag(zt1)).gt.epsocc) then
          !zt1=t1*zt1
          do igp=1,ngknr1
            ifg=igfft(igkignr1(igp))
            wfir1(ifg,ispn,j)=wfir2(ifg,ispn,j)+zt1*evecfv1(igp,ist,1)
          end do
        end if
      end do
    end do
    i=0
    do ispn=1,nspinor
      do ist=1,nstfv
        i=i+1
        !zt1=evecsv2(i,j)	
        if (abs(dble(zt1))+abs(aimag(zt1)).gt.epsocc) then
          zt1=t1*zt1
          do igp=1,ngknr2
            ifg=igfft(igkignr2(igp))
            wfir2(ifg,ispn,j)=wfir2(ifg,ispn,j)+zt1*evecfv2(igp,ist,1)
          end do
        end if
      end do
    end do
  else
! spin-unpolarised wavefunction
    do igp=1,ngknr1
      ifg=igfft(igkignr1(igp))
      wfir1(ifg,1,j)=evecfv1(igp,j,1)
    end do
    do igp=1,ngknr2
      ifg=igfft(igkignr2(igp))
      wfir2(ifg,1,j)=evecfv2(igp,j,1)
    end do
  end if
! Fourier transform wavefunction to real-space
  do ispn=1,nspinor
    call zfftifc(3,ngrid,1,wfir1(:,ispn,j))
    call zfftifc(3,ngrid,1,wfir2(:,ispn,j))
  end do
end do

allocate(zrhofc_tmp(ngvecme,nmemax))
zrhofc_tmp=dcmplx(0.d0,0.d0)
call idxbos(nme0,mpi_dims(2),mpi_x(2)+1,idx0,bs)
i1=idx0+1
i2=idx0+bs
do i=i1,i2
  ist1=ime0(1,i)
  ist2=ime0(2,i)
  ! interstitial part
  if (spinpol) then
  ! spin-polarised
    do ir=1,ngrtot
      zrhoir(ir)=(conjg(wfir1(ir,1,ist1))*wfir2(ir,1,ist2) + &
        conjg(wfir1(ir,2,ist1))*wfir2(ir,2,ist2))*cfunir(ir)
    end do
  else
  ! spin-unpolarised
    do ir=1,ngrtot
      zrhoir(ir)=conjg(wfir1(ir,1,ist1))*wfir2(ir,1,ist2)*cfunir(ir)
    end do
  end if
!  do ir=1,ngrtot
!    zrhoir(ir)=zrhoir(ir)*cfunir(ir)*omega
!  enddo
  call zfftifc(3,ngrid,-1,zrhoir)
  do ig=1,ngvecme
    zrhofc_tmp(ig,i)=zrhofc_tmp(ig,i)+zrhoir(igfft1(ig))
  enddo
enddo

if (mpi_dims(2).gt.1) then
  call d_reduce_cart(comm_cart_010,.false.,zrhofc_tmp,2*ngvecme*nmemax)
endif
zrhofc0=zrhofc0+zrhofc_tmp



deallocate(wfir1,wfir2,zrhoir,zrhofc_tmp)



endif
!
!! interstitial part
!do ir=1,ngrtot
!  zrhoir(ir)=zrhoir(ir)*cfunir(ir)*omega
!enddo
!call zfftifc(3,ngrid,-1,zrhoir)
!do ig=1,ngvec_me
!  zrhofc(ig,2)=zrhoir(igfft1(ig))
!  zrhofc(ig,3)=zrhofc(ig,1)+zrhofc(ig,2)
!enddo
return
end