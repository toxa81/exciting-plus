complex(8) function zfinp_(zfmt1,zfmt2,zfir1,zfir2)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   tsh   : .true. if the muffin-tin functions are in spherical harmonics
!           (in,logical)
!   zfmt1 : first complex function in spherical harmonics/coordinates for all
!           muffin-tins (in,complex(lmmaxvr,nrcmtmax,natmtot))
!   zfmt2 : second complex function in spherical harmonics/coordinates for all
!           muffin-tins (in,complex(lmmaxvr,nrcmtmax,natmtot))
!   zfir1 : first complex interstitial function in real-space
!           (in,complex(ngrtot))
!   zfir2 : second complex interstitial function in real-space
!           (in,complex(ngrtot))
! !DESCRIPTION:
!   Calculates the inner product of two complex fuctions over the entire unit
!   cell. The muffin-tin functions should be stored on the coarse radial grid
!   and have angular momentum cut-off {\tt lmaxvr}. In the intersitial region,
!   the integrand is multiplied with the characteristic function, to remove the
!   contribution from the muffin-tin. See routines {\tt zfmtinp} and
!   {\tt gencfun}.
!
! !REVISION HISTORY:
!   Created July 2004 (Sharma)
!EOP
!BOC
implicit none
! arguments
complex(8), intent(in) :: zfmt1(lmmaxvr,nrmtmax,natmtot)
complex(8), intent(in) :: zfmt2(lmmaxvr,nrmtmax,natmtot)
complex(8), intent(in) :: zfir1(ngrtot)
complex(8), intent(in) :: zfir2(ngrtot)
! local variables
integer is,ia,ias,ir
complex(8) zsumir,zsummt
complex(8) zf1(nrmtmax)
complex(8), external :: zdotc

! interstitial contribution
zsumir=zzero
do ir=1,ngrtot
  zsumir=zsumir+cfunir(ir)*conjg(zfir1(ir))*zfir2(ir)
end do
zsumir=zsumir*omega/dble(ngrtot)
! muffin-tin contribution
zsummt=zzero
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ir=1,nrmt(is)
      zf1(ir)=zdotc(lmmaxvr,zfmt1(1,ir,ias),1,zfmt2(1,ir,ias),1)*spr(ir,is)**2
    enddo
    do ir=1,nrmt(is)-1
      zsummt=zsummt+0.5d0*(zf1(ir)+zf1(ir+1))*(spr(ir+1,is)-spr(ir,is))
    enddo
  end do
end do
zsummt=zsummt*fourpi/dble(lmmaxvr)
zfinp_=zsumir+zsummt
return
end function
