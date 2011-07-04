
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genpmat
! !INTERFACE:
subroutine genpmat(ngp,igpig,vgpc,apwalm,evecfv,evecsv,pmat)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   ngp    : number of G+p-vectors (in,integer)
!   igpig  : index from G+p-vectors to G-vectors (in,integer(ngkmax))
!   vgpc   : G+p-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   evecfv : first-variational eigenvector (in,complex(nmatmax,nstfv))
!   evecsv : second-variational eigenvectors (in,complex(nstsv,nstsv))
!   pmat   : momentum matrix elements (out,complex(3,nstsv,nstsv))
! !DESCRIPTION:
!   Calculates the momentum matrix elements
!   $$ p_{ij}=\langle\Psi_{i,{\bf k}}|-i\nabla|\Psi_{j,{\bf k}}\rangle. $$
!
! !REVISION HISTORY:
!   Created November 2003 (Sharma)
!   Fixed bug found by Juergen Spitaler, September 2006 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
complex(8), intent(out) :: pmat(3,nstsv,nstsv)
! local variables
integer ispn,is,ia,ist,jst
integer i,l,igp,ifg,ir,ias,lm,ic,io,wfsz
complex(8) zt1,zt2,zsum
! allocatable arrays
complex(8), allocatable :: wfmt(:,:)
complex(8), allocatable :: gwfmt(:,:,:)
complex(8), allocatable :: gwfir(:)
complex(8), allocatable :: pm(:,:,:)
complex(8), allocatable :: wffvmt(:,:,:,:)
complex(8), allocatable :: wftmp1(:,:)
complex(8), allocatable :: wftmp2(:,:)
complex(8), allocatable :: zm1(:,:)
complex(8), allocatable :: zm2(:,:)
complex(8), allocatable :: zv1(:,:)

allocate(wffvmt(lmmaxapw,nufrmax,natmtot,nstfv))
call genwffvmt(lmaxapw,lmmaxapw,ngp,evecfv,apwalm,wffvmt)

wfsz=lmmaxapw*nufrmax*natmtot+ngp
allocate(pm(nstfv,nstfv,3))
pm(:,:,:)=0.d0
allocate(wfmt(lmmaxapw,nrmtmax))
allocate(gwfmt(lmmaxapw,nrmtmax,3))
allocate(gwfir(ngrtot))
allocate(wftmp1(wfsz,nstfv))
allocate(wftmp2(wfsz,3))
allocate(zv1(nstfv,3))
! loop over |ket> states 
do ist=1,nstfv
  wftmp1=zzero
  wftmp2=zzero
! muffin-tin part of \grad |ket>
  do ias=1,natmtot
    is=ias2is(ias)
    ia=ias2ia(ias)
    ic=ias2ic(ias)
! calculate the wavefunction
    call wavefmt(1,lmaxapw,is,ia,ngp,apwalm,evecfv(:,ist),lmmaxapw,wfmt)
! calculate the gradient
    call gradzfmt(lmaxapw,nrmt(is),spr(:,is),lmmaxapw,nrmtmax,wfmt,gwfmt)
    do i=1,3
      do lm=1,lmmaxapw
        l=lm2l(lm)
        do io=1,nufr(l,is)
          zsum=zzero
          ir=1
          zt1=gwfmt(lm,ir,i)*ufr(ir,l,io,ic)*(spr(ir,is)**2)
          do ir=2,nrmt(is)
            zt2=gwfmt(lm,ir,i)*ufr(ir,l,io,ic)*(spr(ir,is)**2)
            zsum=zsum+0.5d0*(zt2+zt1)*(spr(ir,is)-spr(ir-1,is))
            zt1=zt2
          enddo
          wftmp2((ias-1)*nufrmax*lmmaxapw+(io-1)*lmmaxapw+lm,i)=zsum
        enddo !io
      enddo !lm
    enddo !i
  enddo !ias
! interstitial part of \grad |ket>
  do i=1,3
    gwfir=zzero
    do igp=1,ngp
      ifg=igfft(igpig(igp))
      zt1=evecfv(igp,ist)
! calculate the gradient, i.e. multiply by i(G+p)
      gwfir(ifg)=vgpc(i,igp)*cmplx(-dimag(zt1),dreal(zt1),8)
    end do !igp
! Fourier transform gradient to real-space, multiply by step function
!  and transform back to G-space
    call zfftifc(3,ngrid,1,gwfir)
    do ir=1,ngrtot
      gwfir(ir)=gwfir(ir)*cfunir(ir)
    end do
    call zfftifc(3,ngrid,-1,gwfir)
    do igp=1,ngp
      wftmp2(natmtot*nufrmax*lmmaxapw+igp,i)=gwfir(igfft(igpig(igp)))
    end do
  end do !i
! collect <bra| states
  do jst=1,ist
! muffin tin part
    call memcopy(wffvmt(1,1,1,jst),wftmp1(1,jst),natmtot*nufrmax*lmmaxapw*16)
! interstitial part
    call memcopy(evecfv(1,jst),wftmp1(natmtot*nufrmax*lmmaxapw+1,jst),ngp*16)
  enddo !jst
  call zgemm('C','N',ist,3,wfsz,zone,wftmp1,wfsz,wftmp2,wfsz,zzero,zv1,nstfv)
  do jst=1,ist
    do i=1,3
      pm(jst,ist,i)=zv1(jst,i)
    enddo
  enddo
end do
deallocate(wfmt,gwfmt,wffvmt,wftmp1,wftmp2,zv1)
deallocate(gwfir)
! multiply by -i and set lower triangular part
do ist=1,nstfv
  do jst=ist,nstfv
    pm(ist,jst,:)=-zi*pm(ist,jst,:)
    pm(jst,ist,:)=conjg(pm(ist,jst,:))
  end do
end do
! compute the second-variational momentum matrix elements
if (tevecsv) then
  pmat=zzero
  allocate(zm1(nstsv,nstfv))
  allocate(zm2(nstsv,nstsv))
  do i=1,3
    do ispn=1,nspinor
      call zgemm('C','N',nstsv,nstfv,nstfv,zone,evecsv((ispn-1)*nstfv+1,1),&
        nstsv,pm(1,1,i),nstfv,zzero,zm1,nstsv)
      call zgemm('N','N',nstsv,nstsv,nstfv,zone,zm1,nstsv,&
        evecsv((ispn-1)*nstfv+1,1),nstsv,zzero,zm2,nstsv)
      pmat(i,:,:)=pmat(i,:,:)+zm2(:,:)
    enddo
  enddo
  deallocate(zm1)
  deallocate(zm2)
else
  do ist=1,nstfv
    do jst=1,nstfv
      pmat(:,ist,jst)=pm(ist,jst,:)
    enddo
  enddo
end if
deallocate(pm)
return
end subroutine
!EOC

