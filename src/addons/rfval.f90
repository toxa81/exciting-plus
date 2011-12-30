subroutine rfval(vrcnr,lmax,ld,rfmt,rfit,val)
use modmain
implicit none
! arguments
real(8), intent(in) :: vrcnr(3)
integer, intent(in) :: lmax
integer, intent(in) :: ld
real(8), intent(in) :: rfmt(ld,nrmtmax,natmtot)
complex(8), intent(in) :: rfit(ngrtot)
real(8), intent(out) :: val
! local variables
integer l,m,lm,ig,ifg,i,j
real(8) tp(2),t1,t2
integer ntr(3),ngloc,igloc
integer ia,is,ias
integer ir0
real(8) vrc0(3),r0
! automatic arrays
real(8) ya(nprad),c(nprad)
real(8) rlm((lmax+1)**2)
logical l2,tmt
real(8) a,b,f0,f1,rmin
integer ir
real(8) val_mt,val_it
! external functions
real(8), external :: polynom
logical, external :: vrinmt
tmt=vrinmt(vrcnr,is,ia,ntr,vrc0,ir0,r0)
val_mt=0.d0
val_it=0.d0
rmin=0.1d0
if (tmt) then
  if (mpi_grid_root()) then
! we will replace the value of the function near the center of the MT sphere
!  with exponent matched by value and first derivative
    l2=.false.
    if (r0.lt.rmin) then
      l2=.true.
! find the first point outside the small sphere around the origin
      do ir=1,nrmt(is)
        if (spr(ir,is).gt.rmin) then
          ir0=ir-nprad/2
          goto 10
        endif
      enddo
 10   continue
    endif
    ias=idxas(ia,is)
    call sphcrd(vrc0,t1,tp)
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
          f0=polynom(0,nprad,spr(ir0,is),ya,c,rmin)
          f1=polynom(1,nprad,spr(ir0,is),ya,c,rmin)
          b=f1/(f0*2*rmin)
          a=f0/exp(b*rmin**2)
          t2=a*exp(b*r0**2)
        else
          t2=polynom(0,nprad,spr(ir0,is),ya,c,r0)
        endif
        val_mt=val_mt+t2*rlm(lm)
      end do
    enddo
  endif
else
  ngloc=mpi_grid_map(ngvec,dim_k)
! otherwise use interstitial function
  do igloc=1,ngloc
    ig=mpi_grid_map(ngvec,dim_k,loc=igloc)
    ifg=igfft(ig)
    t1=vgc(1,ig)*vrcnr(1)+vgc(2,ig)*vrcnr(2)+vgc(3,ig)*vrcnr(3)
    val_it=val_it+dble(rfit(ifg)*cmplx(cos(t1),sin(t1),8))
  end do
endif
val=val_mt+val_it
return
end subroutine
