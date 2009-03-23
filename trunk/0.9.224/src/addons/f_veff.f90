subroutine f_veff(vrc,val)
use modmain
implicit none
real(8), intent(in) :: vrc(3)
real(8), intent(out) :: val

integer ntr(3)
real(8) vr0l(3)

call getntr(avec,vrc,ntr,vr0l)
call rfarray(lmaxvr,lmmaxvr,veffmt,veffir,1,vr0l,val)

return
end

subroutine f_veff_p(vrc,zfft_vir,val)
use modmain
implicit none
real(8), intent(in) :: vrc(3)
complex(8), intent(in) :: zfft_vir(ngrtot) 
real(8), intent(out) :: val

integer ntr(3)
real(8) vr0l(3),v1(3)

call getntr(avec,vrc,ntr,vr0l)
call r3mv(avec,vr0l,v1)
call rfarray_p(lmaxvr,lmmaxvr,veffmt,zfft_vir,v1,val)

return
end

subroutine rfarray_p(lmax,ld,rfmt,zfft,vrc,fp)
use modmain
implicit none
! arguments
integer, intent(in) :: lmax
integer, intent(in) :: ld
real(8), intent(in) :: rfmt(ld,nrmtmax,natmtot)
complex(8), intent(in) :: zfft(ngrtot)
real(8), intent(in) :: vrc(3)
real(8), intent(out) :: fp
! local variables
integer ia,is,ias
integer ir0
integer l,m,lm,ig,ifg,i,j
real(8) tp(2),sum,t1,t2,sum2
integer idx0,bs,ntr(3)
real(8) vr0(3),r0
! automatic arrays
real(8) ya(nprad),c(nprad)
! allocatable arrays
real(8), allocatable :: rlm(:)
logical l1
logical, external :: vrinmt

! external functions
real(8), external :: polynom
allocate(rlm((lmax+1)**2))

l1=vrinmt(vrc,is,ia,ntr,vr0,ir0,r0)
sum=0.d0
sum2=0.d0
if (l1) then
  if (iproc.eq.0) then
    ias=idxas(ia,is)
    call sphcrd(vr0,t1,tp)
    call genrlm(lmax,tp,rlm)
    do l=0,lmax
      do m=-l,l
        lm=idxlm(l,m)
        do j=1,nprad
          i=ir0+j-1
          ya(j)=rfmt(lm,i,ias)
        end do
        t2=polynom(0,nprad,spr(ir0,is),ya,c,r0)
        sum=sum+t2*rlm(lm)
      end do
    enddo
  endif
else
  call idxbos(ngvec,nproc,iproc+1,idx0,bs)
! otherwise use interstitial function
  do ig=idx0+1,idx0+bs
    ifg=igfft(ig)
    t1=vgc(1,ig)*vrc(1)+vgc(2,ig)*vrc(2)+vgc(3,ig)*vrc(3)
    sum2=sum2+dble(zfft(ifg)*cmplx(cos(t1),sin(t1),8))
  end do
  call dsync(sum2,1,.true.,.false.)
endif
if (iproc.eq.0) then
  fp=sum+sum2
endif
call barrier(comm_world)
deallocate(rlm)
return
end subroutine

