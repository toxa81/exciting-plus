subroutine getevecfd(vpl,vgpl,evecfd)
use modmain
implicit none
real(8), intent(in) :: vpl(3)
real(8), intent(in) :: vgpl(3,ngkmax)
complex(8), intent(out) :: evecfd(nspinor*nmatmax,nstsv)
!
integer isym,ik,n,ngkmax_,nmatmax_,nspinor_,nstsv_,lspl,ilspl,j
integer igk,ig,i,offs,is,ia,ias,igp,ilo,jas,ispn,ja,l,lm
integer lspn
real(8) vkl_(3),t1,v(3),si(3,3),det,th
complex(8) zt1,zt2,su2(2,2)
real(8), allocatable :: vgkl_(:,:)
integer, allocatable :: igkig_(:)
complex(8), allocatable :: evectmp(:,:)
!
allocate(vgkl_(3,ngkmax))
allocate(igkig_(ngkmax))
! find the equivalent k-point number and crystal symmetry element
call findkpt(vpl,isym,ik)
! find the record length
inquire(iolength=n) vkl_,vgkl_,igkig_,ngkmax_,nmatmax_,nspinor_,&
  nstsv_,evecfd
open(70,file=trim(scrpath)//'EVECFD'//trim(filext),action='READ', &
  form='UNFORMATTED',access='DIRECT',recl=n)
read(70,rec=ik) vkl_,vgkl_,igkig_,ngkmax_,nmatmax_,nspinor_,&
  nstsv_,evecfd
close(70)
t1=abs(vkl(1,ik)-vkl_(1))+abs(vkl(2,ik)-vkl_(2))+abs(vkl(3,ik)-vkl_(3))
if (t1.gt.epslat) then
  write(*,*)
  write(*,'("Error(getevecfd): differing vectors for k-point ",I8)') ik
  write(*,'(" current    : ",3G18.10)') vkl(:,ik)
  write(*,'(" EVECFD.OUT : ",3G18.10)') vkl_
  write(*,*)
  call pstop
end if
if (ngkmax.ne.ngkmax_) then
  write(*,*)
  write(*,'("Error(getevecfd): differing ngkmax for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') ngkmax
  write(*,'(" EVECFD.OUT : ",I8)') ngkmax_
  write(*,*)
  call pstop
end if
if (nmatmax.ne.nmatmax_) then
  write(*,*)
  write(*,'("Error(getevecfd): differing nmatmax for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nmatmax
  write(*,'(" EVECFD.OUT : ",I8)') nmatmax_
  write(*,*)
  call pstop
end if
if (nstsv.ne.nstsv_) then
  write(*,*)
  write(*,'("Error(getevecsv): differing nstsv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nstsv
  write(*,'(" EVECFV.OUT : ",I8)') nstsv_
  write(*,*)
  call pstop
end if
if (nspinor.ne.nspinor_) then
  write(*,*)
  write(*,'("Error(getevecfv): differing nspnfv for k-point ",I8)') ik
  write(*,'(" current    : ",I8)') nspinor
  write(*,'(" EVECFV.OUT : ",I8)') nspinor_
  write(*,*)
  call pstop
end if
! if p = k then return
t1=abs(vpl(1)-vkl(1,ik))+abs(vpl(2)-vkl(2,ik))+abs(vpl(3)-vkl(3,ik))
if (t1.lt.epslat) then
  deallocate(vgkl_,igkig_)
  return
end if
allocate(evectmp(nspinor*nmatmax,nstsv))   
! index to spatial rotation in lattice point group
lspl=lsplsymc(isym)
! the inverse of the spatial symmetry rotates k into p
ilspl=isymlat(lspl)
si(:,:)=dble(symlat(:,:,ilspl))
!-----------------------------------------------!
!     translate and rotate APW coefficients     !
!-----------------------------------------------!
do ispn=1,nspinor
  offs=(ispn-1)*nmatmax
  do igk=1,ngk(1,ik)
    ig=igkig_(igk)
    v(:)=dble(ivg(:,ig))
    t1=-twopi*dot_product(v(:),vtlsymc(:,isym))
    zt1=cmplx(cos(t1),sin(t1),8)
    evectmp(offs+igk,:)=zt1*evecfd(offs+igk,:)
  end do
! inverse rotation used because transformation is passive
  do igk=1,ngk(1,ik)
    call r3mtv(si,vgkl_(:,igk),v)
    do igp=1,ngk(1,ik)
      t1=abs(v(1)-vgpl(1,igp)) &
        +abs(v(2)-vgpl(2,igp)) &
        +abs(v(3)-vgpl(3,igp))
      if (t1.lt.epslat) then
        evecfd(offs+igp,:)=evectmp(offs+igk,:)
        goto 10
      end if
    end do
10 continue
  end do
end do
!---------------------------------------------------------!
!     translate and rotate local-orbital coefficients     !
!---------------------------------------------------------!
if (nlotot.gt.0) then
! rotate k-point by inverse symmetry matrix
  call r3mtv(si,vkl(:,ik),v)
  do ispn=1,nspinor
    offs=(ispn-1)*nmatmax
! make a copy of the local-orbital coefficients
    do i=ngk(1,ik)+1,nmat(1,ik)
      evectmp(offs+i,:)=evecfd(offs+i,:)
    end do
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
! rotate local orbitals (active transformation)
        do ilo=1,nlorb(is)
          l=lorbl(ilo,is)
          lm=idxlm(l,-l)
          i=ngk(1,ik)+idxlo(lm,ilo,ias)
          j=ngk(1,ik)+idxlo(lm,ilo,jas)
          call rotzflm(symlatc(:,:,lspl),l,l,nstsv,nspinor*nmatmax,&
            evectmp(offs+j,1),evecfd(offs+i,1))
          evecfd(offs+i:offs+i+2*l,:)=zt1*evecfd(offs+i:offs+i+2*l,:)
        end do
      end do
    end do
  end do
end if
! apply spin rotation if required
if (spinpol) then
! index to global spin rotation in lattice point group
  lspn=lspnsymc(isym)
! find the SU(2) representation of the spin rotation matrix
  call rotaxang(epslat,symlatc(:,:,lspn),det,v,th)
  call axangsu2(v,th,su2)
! apply SU(2) matrix to spinor wavefunctions (active transformation)
  do j=1,nstsv
    do i=1,nmat(1,ik)
      zt1=evecfd(i,j)
      zt2=evecfd(i+nmatmax,j)
      evecfd(i,j)=su2(1,1)*zt1+su2(1,2)*zt2
      evecfd(i+nmatmax,j)=su2(2,1)*zt1+su2(2,2)*zt2
    end do
  end do
end if
deallocate(evectmp)
deallocate(vgkl_,igkig_)

return
end subroutine
