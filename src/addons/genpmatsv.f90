subroutine genpmatsv(ngp,igpig,vgpc,wfsvmt,wfsvit,pmat)
use modmain
use mod_util
implicit none
! arguments
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(in) :: wfsvmt(lmmaxapw,nufrmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfsvit(ngkmax,nspinor,nstsv)
complex(8), intent(out) :: pmat(3,nstsv,nstsv)
! local variables
integer ispn,is,ia,ist,jst,n
integer i,l,igp,ifg,ir,ias,lm,ic,io,wfsize,wfsizemax,ig
complex(8) zt1,zt2,zsum
integer idx(3)
! allocatable arrays
complex(8), allocatable :: wfmt(:,:)
complex(8), allocatable :: gwfmt(:,:,:)
complex(8), allocatable :: gwfir(:)
complex(8), allocatable :: wffvmt(:,:,:,:)
complex(8), allocatable :: wftmp1(:,:)
complex(8), allocatable :: wftmp2(:,:)
complex(8), allocatable :: zv1(:,:)
complex(8), allocatable :: zf(:)

wfsizemax=nspinor*(lmmaxapw*nufrmax*natmtot+ngp)
allocate(wfmt(lmmaxapw,nrmtmax))
allocate(gwfmt(lmmaxapw,nrmtmax,3))
allocate(gwfir(ngrtot))
allocate(wftmp1(wfsizemax,nstsv))
allocate(wftmp2(wfsizemax,3))
allocate(zv1(nstsv,3))
allocate(zf(nrmtmax))
! collect <bra| states
wftmp1=zzero
i=0
do ispn=1,nspinor
  do ias=1,natmtot
    is=ias2is(ias)
    do l=0,lmaxapw
      do io=1,nufr(l,is)
        do lm=l**2+1,(l+1)**2
          i=i+1
          wftmp1(i,:)=wfsvmt(lm,io,ias,ispn,:)
        enddo
      enddo
    enddo
  enddo
  do ig=1,ngp
    i=i+1
    wftmp1(i,:)=wfsvit(ig,ispn,:)
  enddo
enddo !ispn
wfsize=i
! loop over |ket> states 
do ist=1,nstsv
  wftmp2=zzero
! muffin-tin part of \grad |ket>
  idx=0
  do ispn=1,nspinor
    do ias=1,natmtot
      wfmt(:,:)=zzero
      is=ias2is(ias)
      ic=ias2ic(ias)
      do l=0,lmaxapw
        do io=1,nufr(l,is)
          do lm=l**2+1,(l+1)**2
            wfmt(lm,:)=wfmt(lm,:)+wfsvmt(lm,io,ias,ispn,ist)*ufr(:,l,io,ic)
          enddo
        enddo
      enddo !l
! calculate the gradient
      call gradzfmt(lmaxapw,nrmt(is),spr(:,is),lmmaxapw,nrmtmax,wfmt,gwfmt)
      ! \grad |psi_j'> = \sum_{lm}gwfmt(lm,ir,i)*Y_{lm}(r)
      ! psi_j> = \sum_{lm,o} A_{lm,o} u_{l,o}(r)*Y_{lm}(r)
      ! <psi_j | \grad |psi_{j'}> = \sum_{lm} \sum_{o} A_{lm,o}^{*} <u_{lo}|gwfmt(lm,ir,i)>
      do i=1,3
        do l=0,lmaxapw
          do io=1,nufr(l,is)
            do lm=l**2+1,(l+1)**2
              do ir=1,nrmt(is)
                zf(ir)=gwfmt(lm,ir,i)*ufr(ir,l,io,ic)
              enddo
              idx(i)=idx(i)+1
              wftmp2(idx(i),i)=zintegrate(nrmt(is),spr(1,is),zf)
            enddo !lm
          enddo !io
        enddo !l
      enddo !i
    enddo !ias
! interstitial part of \grad |ket>
    do i=1,3
      gwfir=zzero
      do igp=1,ngp
        zt1=wfsvit(igp,ispn,ist)
! calculate the gradient, i.e. multiply by i(G+p)
        gwfir(igfft(igpig(igp)))=vgpc(i,igp)*cmplx(-dimag(zt1),dreal(zt1),8)
      enddo !igp
! Fourier transform gradient to real-space, multiply by step function
!  and transform back to G-space
      call zfftifc(3,ngrid,1,gwfir)
      do ir=1,ngrtot
        gwfir(ir)=gwfir(ir)*cfunir(ir)
      end do
      call zfftifc(3,ngrid,-1,gwfir)
      do igp=1,ngp
        idx(i)=idx(i)+1
        wftmp2(idx(i),i)=gwfir(igfft(igpig(igp)))
      end do
    end do !i
  enddo !ispn
  call zgemm('C','N',ist,3,wfsize,zone,wftmp1,wfsizemax,wftmp2,wfsizemax,zzero,zv1,nstsv)
  do jst=1,ist
    do i=1,3
      pmat(i,jst,ist)=zv1(jst,i)
    enddo
  enddo
end do
deallocate(wfmt,gwfmt,wftmp1,wftmp2,zv1,zf)
deallocate(gwfir)
! multiply by -i and set lower triangular part
do ist=1,nstsv
  do jst=ist,nstsv
    pmat(:,ist,jst)=-zi*pmat(:,ist,jst)
    pmat(:,jst,ist)=dconjg(pmat(:,ist,jst))
  end do
end do
return
end subroutine
!EOC

