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
complex(8), intent(in) :: wfsvit1(ngkmax,nspinor,nstsv)
complex(8), intent(in) :: wfsvit2(ngkmax,nspinor,nstsv)
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
real(8) v1(3),v2(3),tp3g(2),len3g
complex(8) sfac3g(natmtot),zt1
complex(8), external :: zdotu,zdotc

allocate(zrhofc_tmp(ngvecme,nmemax))
zrhofc_tmp=dcmplx(0.d0,0.d0)

if (.not.lfftit) then
  allocate(a(ngknr1,nstsv,nspinor))
  allocate(mit(ngknr1,ngknr2))
  
  call idxbos(ngvecme,mpi_dims(2),mpi_x(2)+1,idx0,bs)
  idx_g1=idx0+1
  idx_g2=idx0+bs
  
  do ig=idx_g1,idx_g2
    call genpwit(ngknr1,ngknr2,igkignr1,igkignr2,-(ivg(:,ig+gvecme1-1)+ivg(:,igkq)),mit)
    a=zzero
    do ispn=1,nspinor
      do i=1,nstsv
        call zgemv('N',ngknr1,ngknr2,zone,mit,ngknr1,wfsvit2(1,ispn,i),1,zzero,a(1,i,ispn),1)
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
!        do ig1=1,ngknr1
!          zrhofc_tmp(ig,i)=zrhofc_tmp(ig,i)+dconjg(wfsvit1(ig1,ispn,ist1))*a(ig1,ist2,ispn2)
!        enddo
! this should also work
       zrhofc_tmp(ig,i)=zrhofc_tmp(ig,i)+&
         zdotc(ngknr1,wfsvit1(1,ispn,ist1),1,a(1,ist2,ispn2),1)
      enddo
    enddo
  enddo !ig  
  deallocate(mit,a)
else
  allocate(wfir1(ngrtot,nspinor,nstsv))
  allocate(wfir2(ngrtot,nspinor,nstsv))
  allocate(zrhoir(ngrtot))
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
          zt1=evecsv2(i,j)  
          if (abs(dble(zt1))+abs(aimag(zt1)).gt.epsocc) then
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
  call idxbos(nme0,mpi_dims(2),mpi_x(2)+1,idx0,bs)
  i1=idx0+1
  i2=idx0+bs
  do i=i1,i2
    ist1=ime0(1,i)
    ist2=ime0(2,i)
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
    call zfftifc(3,ngrid,-1,zrhoir)
    do ig=1,ngvecme
      zrhofc_tmp(ig,i)=zrhofc_tmp(ig,i)+zrhoir(igfft1(ig))
    enddo
  enddo
  deallocate(wfir1,wfir2,zrhoir)
endif

if (mpi_dims(2).gt.1) then
  call d_reduce_cart(comm_cart_010,.false.,zrhofc_tmp,2*ngvecme*nmemax)
endif
zrhofc0=zrhofc0+zrhofc_tmp
deallocate(zrhofc_tmp)

return
end