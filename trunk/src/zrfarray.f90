
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zrfarray(lmax,ld,zrfmt,zrfir,np,vpl,zfp)
use modmain
implicit none
! arguments
integer, intent(in) :: lmax
integer, intent(in) :: ld
complex(8), intent(in) :: zrfmt(ld,nrcmtmax,natmtot)
complex(8), intent(in) :: zrfir(ngrtot)
integer, intent(in) :: np
real(8), intent(in) :: vpl(3,np)
complex(8), intent(out) :: zfp(np)
! local variables
integer ia,is,ias,ip,iv(3)
integer i1,i2,i3,ir0,ir,np2
integer l,m,lm,ig,ifg,i,j
real(8) rmt2,r,tp(2),t1
real(8) v1(3),v2(3),v3(3),v4(3),v5(3)
complex(8) :: sum,t2
! automatic arrays
real(8) ya(nprad),c(nprad)
! allocatable arrays
complex(8), allocatable :: ylm(:)
complex(8), allocatable :: zfft(:)
! external functions
real(8) polynom
external polynom
allocate(ylm((lmax+1)**2))
allocate(zfft(ngrtot))
np2=nprad/2
! Fourier transform rfir to G-space
zfft(:)=zrfir(:)
call zfftifc(3,ngrid,-1,zfft)

! begin loop over all points
do ip=1,np
  v2(:)=vpl(:,ip)
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
              call genylm(lmax,tp,ylm)
	      r = max(r,rcmt(1,is))+1.d-10
              do ir=1,nrcmt(is)-1
                if (r.gt.rcmt(ir,is).and.r.le.rcmt(ir+1,is)) then
                  sum=dcmplx(0.d0,0.d0)
		  
                  do l=0,lmax
                    do m=-l,l
                      lm=idxlm(l,m)
		      t2 = zrfmt(lm,ir,ias) + &
		        ((zrfmt(lm,ir+1,ias)-zrfmt(lm,ir,ias))/(rcmt(ir+1,is)-rcmt(ir,is)))* &
		        (r-rcmt(ir,is))
                      sum=sum+t2*ylm(lm)
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
! otherwise use interstitial function
  sum=dcmplx(0.d0,0.d0)
  do ig=1,ngvec
    ifg=igfft(ig)
    t1=vgc(1,ig)*v1(1)+vgc(2,ig)*v1(2)+vgc(3,ig)*v1(3)
    sum=sum+zfft(ifg)*cmplx(cos(t1),sin(t1),8)
  end do
10 continue
  zfp(ip)=sum
end do
deallocate(ylm,zfft)
return
end subroutine
