subroutine f_veff(vrc,val)
use modmain
implicit none
real(8), intent(in) :: vrc(3)
real(8), intent(out) :: val

integer ntr(3)
real(8) vr0l(3),v1(3)
complex(8), allocatable :: zfft(:)
!allocate(zfft(ngrtot))
!zfft(:)=veffir(:)
!call zfftifc(3,ngrid,-1,zfft)
call getntr(avec,vrc,ntr,vr0l)
call r3mv(avec,vr0l,v1)
call rfarray_p(lmaxvr,lmmaxvr,veffmt,veffir_zfft,v1,val)
!deallocate(zfft)
!call getntr(avec,vrc,ntr,vr0l)
!call rfarray(lmaxvr,lmmaxvr,veffmt,veffir,1,vr0l,val)
return
end

subroutine f_veff_p(vrc,vmt,vig,val)
use modmain
implicit none
real(8), intent(in) :: vrc(3)
real(8), intent(in) :: vmt(lmmaxvr,nrmtmax,natmtot)
complex(8), intent(in) :: vig(ngrtot) 
real(8), intent(out) :: val

integer ntr(3)
real(8) vr0l(3),v1(3)

call getntr(avec,vrc,ntr,vr0l)
call r3mv(avec,vr0l,v1)
call rfarray_p(lmaxvr,lmmaxvr,vmt,vig,v1,val)

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
integer idx0,bs,ntr(3),ngloc,igloc
real(8) vr0(3),r0
! automatic arrays
real(8) ya(nprad),c(nprad)
! allocatable arrays
real(8), allocatable :: rlm(:)
logical l1,l2
real(8) a,b,f0,f1
logical, external :: vrinmt
integer ir

! external functions
real(8), external :: polynom
allocate(rlm((lmax+1)**2))

l1=vrinmt(vrc,is,ia,ntr,vr0,ir0,r0)
sum=0.d0
sum2=0.d0
if (l1) then
  if (mpi_grid_root()) then
    l2=.false.
    if (r0.lt.0.1) then
      l2=.true.
      do ir=1,nrmt(is)
        if (spr(ir,is).gt.0.1) then
          ir0=ir-nprad/2
          goto 10
        endif
      enddo
 10   continue
    endif
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
        if (l2.and.l.eq.0) then
          do j=1,nprad
            i=ir0+j-1
            ya(j)=rfmt(lm,i,ias)
          end do
          f0=polynom(0,nprad,spr(ir0,is),ya,c,spr(ir0+nprad/2,is))
          f1=polynom(1,nprad,spr(ir0,is),ya,c,spr(ir0+nprad/2,is))
          b=-f1/(f0*2*spr(ir0+nprad/2,is))
          a=-f0/exp(-b*spr(ir0+nprad/2,is)**2)
          t2=-a*exp(-b*r0**2)
        else
          t2=polynom(0,nprad,spr(ir0,is),ya,c,r0)
        endif
        sum=sum+t2*rlm(lm)
      end do
    enddo
  endif
else
  ngloc=mpi_grid_map(ngvec,dim_k)
! otherwise use interstitial function
  do igloc=1,ngloc
    ig=mpi_grid_map(ngvec,dim_k,loc=igloc)
    ifg=igfft(ig)
    t1=vgc(1,ig)*vrc(1)+vgc(2,ig)*vrc(2)+vgc(3,ig)*vrc(3)
    sum2=sum2+dble(zfft(ifg)*cmplx(cos(t1),sin(t1),8))
  end do
  call mpi_grid_reduce(sum2,dims=(/dim_k/))
endif
if (mpi_grid_root()) then
  fp=sum+sum2
endif
deallocate(rlm)
return
end subroutine
