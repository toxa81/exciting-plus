subroutine f_veff(vrc,val)
use modmain
implicit none
real(8), intent(in) :: vrc(3)
real(8), intent(out) :: val

integer ntr(3)
real(8) vr0l(3)

call getntr(vrc,ntr,vr0l)
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
real(8) vr0l(3)

call getntr(vrc,ntr,vr0l)
call rfarray_p(lmaxvr,lmmaxvr,veffmt,zfft_vir,vr0l,val)

return
end

subroutine rfarray_p(lmax,ld,rfmt,zfft,vpl,fp)
use modmain
implicit none
! arguments
integer, intent(in) :: lmax
integer, intent(in) :: ld
real(8), intent(in) :: rfmt(ld,nrmtmax,natmtot)
complex(8), intent(in) :: zfft(ngrtot)
real(8), intent(in) :: vpl(3)
real(8), intent(out) :: fp
! local variables
integer ia,is,ias,iv(3)
integer i1,i2,i3,ir0,ir,np2
integer l,m,lm,ig,ifg,i,j
real(8) rmt2,r,tp(2),sum,t1,t2,sum2
real(8) v1(3),v2(3),v3(3),v4(3),v5(3)
integer idx0,bs
! automatic arrays
real(8) ya(nprad),c(nprad)
! allocatable arrays
real(8), allocatable :: rlm(:)
logical l1
! external functions
real(8) polynom
external polynom
allocate(rlm((lmax+1)**2))
np2=nprad/2
l1=.false.
! begin loop over all points
if (iproc.eq.0) then
  v2(:)=vpl(:)
  call r3frac(epslat,v2,iv)
! convert point to Cartesian coordinates
  call r3mv(avec,v2,v1)
! check if point is in a muffin-tin
  do is=1,nspecies
    rmt2=rmt(is)**2
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      v2(:)=v1(:)-atposc(:,ia,is)
      do i1=-1,1
        v3(:)=v2(:)+dble(i1)*avec(:,1)
        do i2=-1,1
          v4(:)=v3(:)+dble(i2)*avec(:,2)
          do i3=-1,1
            v5(:)=v4(:)+dble(i3)*avec(:,3)
            t1=v5(1)**2+v5(2)**2+v5(3)**2
            if (t1.lt.rmt2) then
              call sphcrd(v5,r,tp)
              call genrlm(lmax,tp,rlm)
              do ir=1,nrmt(is)
                if (spr(ir,is).ge.r) then
                  if (ir.le.np2) then
                    ir0=1
                  else if (ir.gt.nrmt(is)-np2) then
                    ir0=nrmt(is)-nprad+1
                  else
                    ir0=ir-np2
                  end if
                  r=max(r,spr(1,is))
                  sum=0.d0
                  do l=0,lmax
                    do m=-l,l
                      lm=idxlm(l,m)
                      do j=1,nprad
                        i=ir0+j-1
                        ya(j)=rfmt(lm,i,ias)
                      end do
                      t2=polynom(0,nprad,spr(ir0,is),ya,c,r)
                      sum=sum+t2*rlm(lm)
                    end do
                  end do
                  goto 10
                end if
              end do
            end if
          end do
        end do
      end do
    end do
  end do
10 continue
endif
call lsync(l1,1,.true.)
if (.not.l1) then
  call idxbos(ngvec,nproc,iproc+1,idx0,bs)
! otherwise use interstitial function
  sum2=0.d0
  do ig=idx0+1,idx0+bs
    ifg=igfft(ig)
    t1=vgc(1,ig)*v1(1)+vgc(2,ig)*v1(2)+vgc(3,ig)*v1(3)
    sum2=sum2+dble(zfft(ifg)*cmplx(cos(t1),sin(t1),8))
  end do
  call dsync(sum2,1,.true.,.false.)
endif
if (iproc.eq.0) then
  fp=sum+sum2
endif
deallocate(rlm)
return
end subroutine

