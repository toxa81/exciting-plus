! 
! local function (lf) algebra
!
module mod_lf

! maximum lattice translation
integer trmax
! total number of translations
integer ntr
integer ntrloc
! translation vectors in lattice coordinates
integer, allocatable :: vtl(:,:)
! vector -> index map
integer, allocatable :: ivtit(:,:,:)

integer dim_t

contains

subroutine lf_init(trmax_,dim_t_)
use modmain
implicit none
integer, intent(in) :: trmax_ 
integer, intent(in) :: dim_t_
integer i1,i2,i3,n
trmax=trmax_
ntr=(2*trmax+1)**3
if (allocated(vtl)) deallocate(vtl)
allocate(vtl(3,ntr))
if (allocated(ivtit)) deallocate(ivtit)
allocate(ivtit(-trmax:trmax,-trmax:trmax,-trmax:trmax))
n=0
do i1=-trmax,trmax
  do i2=-trmax,trmax
    do i3=-trmax,trmax
      n=n+1
      vtl(:,n)=(/i1,i2,i3/)
      ivtit(i1,i2,i3)=n
    enddo
  enddo
enddo
dim_t=dim_t_
ntrloc=mpi_grid_map(ntr,dim_t)
return
end subroutine

complex(8) function lf_dotlf(tsh,t,f1mt,f1ir,f2mt,f2ir)
use modmain
implicit none
logical, intent(in) :: tsh
integer, intent(in) :: t(3)
complex(8), intent(in) :: f1mt(lmmaxvr,nrmtmax,natmtot,*)
complex(8), intent(in) :: f1ir(ngrtot,*)
complex(8), intent(in) :: f2mt(lmmaxvr,nrmtmax,natmtot,*)
complex(8), intent(in) :: f2ir(ngrtot,*)
complex(8), allocatable :: f2mt_tmp(:,:,:)
complex(8), allocatable :: f2ir_tmp(:)

complex(8) zprod
integer ispn,it,jt,v1(3),v2(3),j,ntstep,ntrloc,itstep,ntloc1
integer jtloc,i,tag
logical l1
complex(8), external :: zfinp_

! compute <f1_0|f2_T>=\int_{-\inf}^{\inf} f1^{*}(r)f2(r-T)dr = 
!   = \sum_{R} \int_{\Omega} f1^{*}(r+R)f2(r+R-T)dr

zprod=zzero
allocate(f2mt_tmp(lmmaxvr,nrmtmax,natmtot))
allocate(f2ir_tmp(ngrtot))

j=0
ntstep=mpi_grid_map(ntr,dim_t,x=j)
ntrloc=mpi_grid_map(ntr,dim_t)
do itstep=1,ntstep
  f2mt_tmp=zzero
  f2ir_tmp=zzero
  do i=0,mpi_grid_size(dim_t)-1
    ntloc1=mpi_grid_map(ntr,dim_t,x=i)
    if (itstep.le.ntloc1) then
      it=mpi_grid_map(ntr,dim_t,x=i,loc=itstep)
      v1(:)=vtl(:,it)
      v2(:)=v1(:)-t(:)
      l1=.false.
      if (v2(1).ge.-trmax.and.v2(1).le.trmax.and.&
          v2(2).ge.-trmax.and.v2(2).le.trmax.and.&
          v2(3).ge.-trmax.and.v2(3).le.trmax) then
        jt=ivtit(v2(1),v2(2),v2(3))
        l1=.true.
        jtloc=mpi_grid_map(ntr,dim_t,glob=jt,x=j)
      endif
      if (l1.and.mpi_grid_x(dim_t).eq.j.and.mpi_grid_x(dim_t).ne.i) then
        tag=(itstep*mpi_grid_size(dim_t)+i)*10
        call mpi_grid_send(f2mt(1,1,1,jtloc),lmmaxvr*nrmtmax*natmtot,&
          (/dim_t/),(/i/),tag)
        call mpi_grid_send(f2ir(1,jtloc),ngrtot,(/dim_t/),(/i/),tag+1)
      endif
      if (l1.and.mpi_grid_x(dim_t).eq.i) then
        if (j.ne.i) then
          tag=(itstep*mpi_grid_size(dim_t)+i)*10
          call mpi_grid_recieve(f2mt_tmp(1,1,1),lmmaxvr*nrmtmax*natmtot,&
            (/dim_t/),(/j/),tag)
          call mpi_grid_recieve(f2ir_tmp(1),ngrtot,(/dim_t/),(/j/),tag+1)
        else
          f2mt_tmp(:,:,:)=f2mt(:,:,:,jtloc)
          f2ir_tmp(:)=f2ir(:,jtloc)
        endif
      endif
    endif
  enddo !
  if (itstep.le.ntrloc) then
    zprod=zprod+zfinp_(tsh,f1mt(1,1,1,itstep),f2mt_tmp,f1ir(1,itstep),&
      f2ir_tmp)
  endif
  call mpi_grid_barrier(dims=(/dim_t/))
enddo
deallocate(f2mt_tmp,f2ir_tmp)
call mpi_grid_reduce(zprod,dims=(/dim_t/))
lf_dotlf=zprod
return
end function

! local function dot bloch function
complex(8) function lf_dotblh(tsh,vpc,f1mt,f1ir,f2mt,f2ir)
use modmain
implicit none
logical, intent(in) :: tsh
real(8), intent(in) :: vpc(3)
complex(8), intent(in) :: f1mt(lmmaxvr,nrmtmax,natmtot,*)
complex(8), intent(in) :: f1ir(ngrtot,*)
complex(8), intent(in) :: f2mt(lmmaxvr,nrmtmax,natmtot)
complex(8), intent(in) :: f2ir(ngrtot)
complex(8) zprod
integer it,itloc
real(8) vtc(3)
complex(8), external :: zfinp_
! <f|psi> = \int dr f^{*}(r) psi(r) =
!   = \sum_R \int_{Omega} dr f^{*}(r+R) psi(r+R) = 
!     \sum_R e^{ikR} \int_{Omega} dr f^{*}(r+R) psi(r)
zprod=zzero
do itloc=1,ntrloc
  it=mpi_grid_map(ntr,dim_t,loc=itloc)
  vtc(:)=vtl(1,it)*avec(:,1)+vtl(2,it)*avec(:,2)+vtl(3,it)*avec(:,3)
  zprod=zprod+exp(zi*dot_product(vpc,vtc))*&
    zfinp_(tsh,f1mt(1,1,1,itloc),f2mt,f1ir(1,itloc),f2ir)
enddo
call mpi_grid_reduce(zprod,dims=(/dim_t/))
lf_dotblh=zprod
return
end function


!subroutine lf_prod(f1mt,f1ir,f2mt,f2ir,f3mt,f3ir)
!use modmain
!implicit none
!complex(8), intent(in) :: f1mt(lmmaxvr,nrmtmax,natmtot,*)
!complex(8), intent(in) :: f1ir(ngrtot,*)
!complex(8), intent(in) :: f2mt(lmmaxvr,nrmtmax,natmtot,*)
!complex(8), intent(in) :: f2ir(ngrtot,*)
!complex(8), intent(out) :: f3mt(lmmaxvr,nrmtmax,natmtot,*)
!complex(8), intent(out) :: f3ir(ngrtot,*)
!
!complex(8), allocatable :: ft1(:,:,:)
!complex(8), allocatable :: ft2(:,:,:)
!integer itrloc,ntrloc
!
!ntrloc=mpi_grid_map(ntr,dim_t)
!allocate(ft1(lmmaxvr,nrmtmax,natmtot))
!allocate(ft2(lmmaxvr,nrmtmax,natmtot))
!do itrloc=1,ntrloc
!  call lf_sht('B',f1mt(1,1,1,itrloc),ft1)
!  call lf_sht('B',f2mt(1,1,1,itrloc),ft2)
!  f3mt(:,:,:,itrloc)=dconjg(ft1(:,:,:))*ft2(:,:,:)
!  f3ir(:,itrloc)=dconjg(f1ir(:,itrloc))*f2ir(:,itrloc)
!enddo
!deallocate(ft1,ft2)
!return 
!end subroutine

! compute f3(r)=alpha*f1(r)*f2(r)+beta*f3(r)
subroutine lf_prod(alpha,f1mt,f1ir,f2mt,f2ir,beta,f3mt,f3ir)
use modmain
implicit none
complex(8), intent(in) :: alpha
complex(8), intent(in) :: f1mt(lmmaxvr,nrmtmax,natmtot,*)
complex(8), intent(in) :: f1ir(ngrtot,*)
complex(8), intent(in) :: f2mt(lmmaxvr,nrmtmax,natmtot,*)
complex(8), intent(in) :: f2ir(ngrtot,*)
complex(8), intent(in) :: beta
complex(8), intent(inout) :: f3mt(lmmaxvr,nrmtmax,natmtot,*)
complex(8), intent(inout) :: f3ir(ngrtot,*)
integer lm1,lm2,lm3,ir,itloc,ias
complex(8), allocatable :: f1mt_(:,:),f2mt_(:,:),f3mt_(:,:)
real(8), external :: gaunt
real(8) t1

allocate(f1mt_(nrmtmax,lmmaxvr))
allocate(f2mt_(nrmtmax,lmmaxvr))
allocate(f3mt_(nrmtmax,lmmaxvr))

do itloc=1,ntrloc
  do ias=1,natmtot
    do lm1=1,lmmaxvr
      f1mt_(:,lm1)=f1mt(lm1,:,ias,itloc)
      f2mt_(:,lm1)=f2mt(lm1,:,ias,itloc)
      f3mt_(:,lm1)=beta*f3mt(lm1,:,ias,itloc)
    enddo
    do lm1=1,lmmaxvr
      do lm2=1,lmmaxvr
        do lm3=1,lmmaxvr
          t1=gaunt(lm2l(lm3),lm2l(lm1),lm2l(lm2),&
                   lm2m(lm3),lm2m(lm1),lm2m(lm2))
          if (abs(t1).gt.1d-8) then
            do ir=1,nrmt(ias2is(ias))
              f3mt_(ir,lm3)=f3mt_(ir,lm3)+alpha*f1mt_(ir,lm1)*f2mt_(ir,lm2)*t1
            enddo
          endif
        enddo
      enddo
    enddo
    do lm3=1,lmmaxvr
      f3mt(lm3,:,ias,itloc)=f3mt_(:,lm3)
    enddo
  enddo !ias
  f3ir(:,itloc)=alpha*f1ir(:,itloc)*f2ir(:,itloc)+beta*f3ir(:,itloc)
enddo !itloc
deallocate(f1mt_,f2mt_,f3mt_)
return
end subroutine




subroutine lf_sht(sht,fmt_in,fmt_out)
use modmain
implicit none
character, intent(in) :: sht
complex(8), intent(in) :: fmt_in(lmmaxvr,nrmtmax,natmtot)
complex(8), intent(out) :: fmt_out(lmmaxvr,nrmtmax,natmtot)
integer ias,is
complex(8), allocatable :: f1(:,:)

allocate(f1(lmmaxvr,nrmtmax))
! backward transform (spherical harmonics to coordinates)
if (sht.eq.'b'.or.sht.eq.'B') then
  do ias=1,natmtot
    is=ias2is(ias)
    f1=fmt_in(:,:,ias)
    call zgemm('N','N',lmmaxvr,nrmt(is),lmmaxvr,zone,zbshtvr, &
      lmmaxvr,f1,lmmaxvr,zzero,fmt_out(:,:,ias),lmmaxvr)
  enddo !ias
endif
! forward transform (spherical coordinates to harmonics)
if (sht.eq.'f'.or.sht.eq.'F') then
  do ias=1,natmtot
    is=ias2is(ias)
    f1=fmt_in(:,:,ias)
    call zgemm('N','N',lmmaxvr,nrmt(is),lmmaxvr,zone,zfshtvr, &
      lmmaxvr,f1,lmmaxvr,zzero,fmt_out(:,:,ias),lmmaxvr)
  enddo !ias
endif
deallocate(f1)
return
end subroutine

!subroutine lf_write(fname,fmt,fir)
!use modmain
!implicit none
!character*(*), intent(in) :: fname
!complex(8), intent(in) :: fmt(lmmaxvr,nrmtmax,natmtot,*)
!complex(8), intent(in) :: fir(ngrtot,*)
!
!integer ntrloc
!
!ntrloc=mpi_grid_map(ntr,dim_t)
!
!end subroutine


end module