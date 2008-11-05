! Copyright (C) 2007 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine rotevecfv(vpl,isym,evecfv,evecfvrot)
use modmain
implicit none
! arguments
real(8), intent(in) :: vpl(3)
integer, intent(in) :: isym
!real(8), intent(in) :: vgpl(3,ngkmax)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(out) :: evecfvrot(nmatmax,nstfv)
! local variables
integer lspl,ilspl
integer ilo,l,m,lm,i
integer ik,igp,igk,ist
integer is,ia,ja,ias,jas
integer recl,nmatmax_,nstfv_,nspnfv_
real(8) vkl_(3),v(3),t1
real(8) si(3,3),sc(3,3)
complex(8) zt1
! allocatable arrays
complex(8), allocatable :: evecfvt(:,:)
complex(8), allocatable :: zflm1(:,:),zflm2(:,:)
! external functions
real(8) r3taxi,r3dot
external r3taxi,r3dot
real(8) vplrot(3),vpcrot(3),v1(3),v2(3)
integer ngprot,ig
integer, allocatable :: igpigrot(:)
real(8), allocatable :: vgplrot(:,:)
real(8), allocatable :: vgpcrot(:,:)
real(8), allocatable :: gpcrot(:)
real(8), allocatable :: tpgpcrot(:,:)
integer, external :: ikloc

!find index of irreducible k-point
call findkpt(vpl,i,ik)
if (i.ne.1) then
  write(*,*)
  write(*,'("Error(rotevecfv) : symmetry operation for irreducible k-point != 1")')
  write(*,*)
  call pstop
endif

lspl=lsplsymc(isym)
ilspl=isymlat(lspl)


!vlprot=\alpha^{-1}(isym)*vpl
si(:,:)=symlat(:,:,ilspl)
call r3frac(epslat,vpl,v1)
call r3mtv(si,v1,v2)
call r3frac(epslat,v2,vplrot)

call r3mv(bvec,vplrot,vpcrot)
allocate(igpigrot(ngkmax))
allocate(vgplrot(3,ngkmax))
allocate(vgpcrot(3,ngkmax))
allocate(gpcrot(ngkmax))
allocate(tpgpcrot(2,ngkmax))

call gengpvec(vplrot,vpcrot,ngprot,igpigrot,vgplrot,vgpcrot,gpcrot,tpgpcrot)

allocate(evecfvt(nmatmax,nstfv))
evecfvt(:,:)=evecfv(:,:)

!-----------------------------------------------!
!     translate and rotate APW coefficients     !
!-----------------------------------------------!
do igk=1,ngk(ik,1)
  ig=igkig(igk,ikloc(ik),1)
  v(:)=dble(ivg(:,ig))
  t1=-twopi*dot_product(v(:),vtlsymc(:,isym))
  zt1=cmplx(cos(t1),sin(t1),8)
  evecfvt(igk,:)=zt1*evecfv(igk,:)
end do
do igk=1,ngk(ik,1)
  call r3mtv(si,vgkl(:,igk,ikloc(ik),1),v)
  do igp=1,ngk(ik,1)
    t1=sum(abs(v(:)-vgplrot(:,igp)))
    if (t1.lt.epslat) then
      evecfvrot(igp,:)=evecfvt(igk,:)
      goto 10
    end if
  end do
10 continue
end do
! return if there are no local-orbitals
if (nlotot.le.0) goto 20
!---------------------------------------------------------!
!     translate and rotate local-orbital coefficients     !
!---------------------------------------------------------!
allocate(zflm1(lolmmax,nstfv),zflm2(lolmmax,nstfv))
! spatial rotation symmetry matrix in Cartesian coordinates
sc(:,:)=symlatc(:,:,lspl)
call r3mtv(si,vkl(:,ik),v)
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
! equivalent atom for this symmetry
    ja=ieqatom(ia,is,isym)
    jas=idxas(ja,is)
! phase factor from translation
    t1=-twopi*dot_product(vkl(:,ik),atposl(:,ja,is))
    zt1=cmplx(cos(t1),sin(t1),8)
    t1=twopi*dot_product(v(:),atposl(:,ia,is))
    zt1=zt1*cmplx(cos(t1),sin(t1),8)
! rotate local orbitals
    do ilo=1,nlorb(is)
      l=lorbl(ilo,is)
      zflm1(:,:)=0.d0
      do m=-l,l
        lm=idxlm(l,m)
        i=ngk(ik,1)+idxlo(lm,ilo,jas)
        zflm1(lm,:)=evecfv(i,:)
      end do
      call rotzflm(sc,l,nstfv,lolmmax,zflm1,zflm2)
      do m=-l,l
        lm=idxlm(l,m)
        i=ngk(ik,1)+idxlo(lm,ilo,ias)
        evecfvrot(i,:)=zt1*zflm2(lm,:)
      end do
    end do
  end do
end do
deallocate(zflm1,zflm2)
20 continue
deallocate(evecfvt,igpigrot,vgplrot,vgpcrot,gpcrot,tpgpcrot)
return
end subroutine

