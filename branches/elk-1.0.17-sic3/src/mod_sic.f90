module mod_sic

! total number of translations
integer ntr
! maximum number of translation vectors
integer, parameter :: maxvtl=1000
! translation vectors in lattice coordinates
integer, allocatable :: vtl(:,:)
! translation vectors in Cartesian coordinates
real(8), allocatable :: vtc(:,:)
! vector -> index map
integer, allocatable :: ivtit(:,:,:)
! translation limits along each lattice vector
integer tlim(2,3)

integer :: ngrloc
integer :: groffs
integer :: nmtloc
integer :: mtoffs
! weights for radial integration
real(8), allocatable :: rmtwt(:)

contains

subroutine sic_copy_mt(tfrwrd,ld,fmt1,fmt2)
use modmain
implicit none
logical, intent(in) :: tfrwrd
integer, intent(in) :: ld
complex(8), intent(inout) :: fmt1(ld,nrmtmax*natmtot)
complex(8), intent(inout) :: fmt2(ld,nmtloc)
if (tfrwrd) then
  fmt2(1:ld,1:nmtloc)=fmt1(1:ld,mtoffs+1:mtoffs+nmtloc)
else
  fmt1(1:ld,mtoffs+1:mtoffs+nmtloc)=fmt2(1:ld,1:nmtloc)
endif
return
end subroutine

subroutine sic_copy_mt_d(tfrwrd,ld,fmt1,fmt2)
use modmain
implicit none
logical, intent(in) :: tfrwrd
integer, intent(in) :: ld
real(8), intent(inout) :: fmt1(ld,nrmtmax*natmtot)
real(8), intent(inout) :: fmt2(ld,nmtloc)
if (tfrwrd) then
  fmt2(1:ld,1:nmtloc)=fmt1(1:ld,mtoffs+1:mtoffs+nmtloc)
else
  fmt1(1:ld,mtoffs+1:mtoffs+nmtloc)=fmt2(1:ld,1:nmtloc)
endif
return
end subroutine

subroutine sic_copy_ir(tfrwrd,fir1,fir2)
use modmain
implicit none
logical, intent(in) :: tfrwrd
complex(8), intent(inout) :: fir1(ngrtot)
complex(8), intent(inout) :: fir2(ngrloc)
if (tfrwrd) then
  fir2(1:ngrloc)=fir1(groffs+1:groffs+ngrloc)
else
  fir1(groffs+1:groffs+ngrloc)=fir2(1:ngrloc)
endif
return
end subroutine

! compute <f1_0|f2_T>=\int_{-\inf}^{\inf} f1^{*}(r)f2(r-T)dr = 
!   = \sum_{R} \int_{\Omega} f1^{*}(r+R)f2(r+R-T)dr
complex(8) function sic_dot_ll(fmt1,fir1,fmt2,fir2,t)
use modmain
implicit none
! arguments
complex(8), intent(in) :: fmt1(lmmaxvr,nmtloc,ntr)
complex(8), intent(in) :: fir1(ngrloc,ntr)
complex(8), intent(in) :: fmt2(lmmaxvr,nmtloc,ntr)
complex(8), intent(in) :: fir2(ngrloc,ntr)
integer, intent(in) :: t(3)
! local variables
integer v1(3),v2(3),jt,ir,is,ias,it,i
complex(8) zdotmt,zdotir,zt1
complex(8), external :: zdotc
zdotmt=zzero
zdotir=zzero
do it=1,ntr
  v1(:)=vtl(:,it)
  v2(:)=v1(:)-t(:)
  if (v2(1).ge.tlim(1,1).and.v2(1).le.tlim(2,1).and.&
      v2(2).ge.tlim(1,2).and.v2(2).le.tlim(2,2).and.&
      v2(3).ge.tlim(1,3).and.v2(3).le.tlim(2,3)) then
    jt=ivtit(v2(1),v2(2),v2(3))
    if (jt.ne.-1) then
      do i=1,nmtloc
        zdotmt=zdotmt+rmtwt(i)*zdotc(lmmaxvr,fmt1(:,i,it),1,&
          fmt2(:,i,jt),1)
      enddo
      do ir=1,ngrloc
        zdotir=zdotir+cfunir(ir+groffs)*dconjg(fir1(ir,it))*fir2(ir,jt)
      enddo
    endif
  endif
enddo
sic_dot_ll=zdotmt+zdotir*omega/dble(ngrtot) 
return
end function

end module