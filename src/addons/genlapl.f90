subroutine genlapl(ngp,igpig,vgpc,wfsvmt,wfsvit,laplsv)
use modmain
implicit none
! arguments
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(in) :: wfsvmt(lmmaxvr,nufrmax,natmtot,nspinor,nstsv)
complex(8), intent(in) :: wfsvit(ngkmax,nspinor,nstsv)
real(8), intent(out) :: laplsv(nstsv)
! local variables
integer ispn,is,ia,j
integer i,l,igp,ifg,ir,ias,lm,ic,io
complex(8) zt1
! allocatable arrays
complex(8), allocatable :: wfmt(:,:)
complex(8), allocatable :: gwfmt(:,:,:)
complex(8), allocatable :: gwfir(:)
complex(8) zf(nrmtmax)
complex(8), external :: zdotc

laplsv=0.d0
allocate(wfmt(lmmaxvr,nrmtmax))
allocate(gwfmt(lmmaxvr,nrmtmax,3))
allocate(gwfir(ngrtot))
do j=1,nstsv
  do ispn=1,nspinor
! muffin tin contribution
    do ias=1,natmtot
      is=ias2is(ias)
      ia=ias2ia(ias)
      ic=ias2ic(ias)
      wfmt=zzero
      do ir=1,nrmt(is)
        do lm=1,lmmaxvr
          l=lm2l(lm)
          do io=1,nufr(l,is)
            wfmt(lm,ir)=wfmt(lm,ir)+wfsvmt(lm,io,ias,ispn,j)*ufr(ir,l,io,ic)
          enddo !io
        enddo !lm
      enddo !ir
! calculate the gradient
      call gradzfmt(lmaxvr,nrmt(is),spr(:,is),lmmaxvr,nrmtmax,wfmt,gwfmt)
! calculate \grad\psi * \grad\psi
      do i=1,3
        do ir=1,nrmt(is)
          zf(ir)=zdotc(lmmaxvr,gwfmt(1,ir,i),1,gwfmt(1,ir,i),1)*(spr(ir,is)**2)
        enddo
        zt1=zzero
        do ir=1,nrmt(is)-1
          zt1=zt1+0.5d0*(zf(ir)+zf(ir+1))*(spr(ir+1,is)-spr(ir,is))
        enddo
        laplsv(j)=laplsv(j)+dreal(zt1)
      enddo !i
    enddo !ias
! interstitial contribution
    do i=1,3
      gwfir=zzero
      do igp=1,ngp
        ifg=igfft(igpig(igp))
        zt1=wfsvit(igp,ispn,j)
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
      zt1=zzero
      do igp=1,ngp
        zt1=zt1+dconjg(gwfir(igfft(igpig(igp))))*gwfir(igfft(igpig(igp)))
      enddo
      laplsv(j)=laplsv(j)+dreal(zt1)
    enddo !i
  enddo !ispn
enddo !j
deallocate(wfmt)
deallocate(gwfmt)
deallocate(gwfir)
return
end subroutine
!EOC

