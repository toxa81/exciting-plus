! 
! local function (lf) algebra
!
module mod_lf

! total number of translations
integer ntr
integer ntrloc
! maximum number of translation vectors
integer, parameter :: maxvtl=1000
! translation vectors in lattice coordinates
integer, allocatable :: vtl(:,:)
! vector -> index map
integer, allocatable :: ivtit(:,:,:)
! translation limits along each lattice vector
integer tlim(2,3)

integer dim_t

contains

complex(8) function lf_dot_lf(tsh,fmt1,fir1,t,fmt2,fir2,tfmt1,tfmt2)
use modmain
implicit none
! arguments
logical, intent(in) :: tsh
integer, intent(in) :: t(3)
logical, intent(in) :: tfmt1(natmtot,ntr)
logical, intent(in) :: tfmt2(natmtot,ntr)
complex(8), intent(in) :: fmt1(lmmaxvr,nrmtmax,natmtot,ntrloc)
complex(8), intent(in) :: fir1(ngrtot,ntrloc)
complex(8), intent(in) :: fmt2(lmmaxvr,nrmtmax,natmtot,ntrloc)
complex(8), intent(in) :: fir2(ngrtot,ntrloc)
! local variables
complex(8), allocatable :: fmt2_(:,:,:)
complex(8), allocatable :: fir2_(:)
!complex(8), external :: zfmtinp_
complex(8) zprod,zt1
integer ispn,it,jt,v1(3),v2(3),j,ntstep,itstep,ntloc1
integer jtloc,i,tag,is,ir,ias
logical l1,l2
complex(8), external :: zfmtinp_
logical, allocatable :: tfmt2_(:)
! compute <f1_0|f2_T>=\int_{-\inf}^{\inf} f1^{*}(r)f2(r-T)dr = 
!   = \sum_{R} \int_{\Omega} f1^{*}(r+R)f2(r+R-T)dr
zprod=zzero
allocate(fmt2_(lmmaxvr,nrmtmax,natmtot))
allocate(fir2_(ngrtot))
allocate(tfmt2_(natmtot))
j=0
ntstep=mpi_grid_map(ntr,dim_t,x=j)
do itstep=1,ntstep
  fmt2_=zzero
  fir2_=zzero
  tfmt2_=.false.
  do i=0,mpi_grid_size(dim_t)-1
    ntloc1=mpi_grid_map(ntr,dim_t,x=i)
    if (itstep.le.ntloc1) then
      it=mpi_grid_map(ntr,dim_t,x=i,loc=itstep)
      v1(:)=vtl(:,it)
      v2(:)=v1(:)-t(:)
      l1=.false.
      if (v2(1).ge.tlim(1,1).and.v2(1).le.tlim(2,1).and.&
          v2(2).ge.tlim(1,2).and.v2(2).le.tlim(2,2).and.&
          v2(3).ge.tlim(1,3).and.v2(3).le.tlim(2,3)) then
        jt=ivtit(v2(1),v2(2),v2(3))
        if (jt.ne.-1) then
          l1=.true.
          jtloc=mpi_grid_map(ntr,dim_t,glob=jt,x=j)
        endif
      endif
      if (l1.and.mpi_grid_x(dim_t).eq.j.and.mpi_grid_x(dim_t).ne.i) then
        tag=(itstep*mpi_grid_size(dim_t)+i)*10
        call mpi_grid_send(fmt2(1,1,1,jtloc),lmmaxvr*nrmtmax*natmtot,&
          (/dim_t/),(/i/),tag)
        call mpi_grid_send(fir2(1,jtloc),ngrtot,(/dim_t/),(/i/),tag+1)
        call mpi_grid_send(tfmt2(1,jt),natmtot,(/dim_t/),(/i/),tag+2)
      endif
      if (l1.and.mpi_grid_x(dim_t).eq.i) then
        if (j.ne.i) then
          tag=(itstep*mpi_grid_size(dim_t)+i)*10
          call mpi_grid_recieve(fmt2_(1,1,1),lmmaxvr*nrmtmax*natmtot,&
            (/dim_t/),(/j/),tag)
          call mpi_grid_recieve(fir2_(1),ngrtot,(/dim_t/),(/j/),tag+1)
          call mpi_grid_recieve(tfmt2_(1),natmtot,(/dim_t/),(/j/),tag+2)
        else
          fmt2_(:,:,:)=fmt2(:,:,:,jtloc)
          fir2_(:)=fir2(:,jtloc)
          tfmt2_(:)=tfmt2(:,jt)
        endif
      endif
    endif !itstep.le.ntloc1
  enddo !i
  if (itstep.le.ntrloc) then
    it=mpi_grid_map(ntr,dim_t,loc=itstep)
    zt1=zzero
    do ir=1,ngrtot
      zt1=zt1+cfunir(ir)*dconjg(fir1(ir,itstep))*fir2_(ir)
    enddo
    zt1=zt1*omega/dble(ngrtot)
    do ias=1,natmtot
      is=ias2is(ias)
      if (tfmt1(ias,it).and.tfmt2_(ias)) then
        zt1=zt1+zfmtinp_(tsh,lmaxvr,nrmt(is),spr(:,is),lmmaxvr,&
          fmt1(1,1,ias,itstep),fmt2_(1,1,ias))
      endif
    enddo
    zprod=zprod+zt1
  endif
  call mpi_grid_barrier(dims=(/dim_t/))
enddo
deallocate(fmt2_,fir2_,tfmt2_)
call mpi_grid_reduce(zprod,dims=(/dim_t/))
lf_dot_lf=zprod
return
end function


! local function dot bloch function
complex(8) function lf_dot_blh(tsh,vpc,f1mt,f1ir,f2mt,f2ir)
use modmain
implicit none
logical, intent(in) :: tsh
real(8), intent(in) :: vpc(3)
complex(8), intent(in) :: f1mt(lmmaxvr,nrmtmax,natmtot,ntrloc)
complex(8), intent(in) :: f1ir(ngrtot,ntrloc)
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
call mpi_grid_reduce(zprod,dims=(/dim_t/),all=.true.)
lf_dot_blh=zprod
return
end function


complex(8) function intgr_zdz(f1mt,f1ir,f2mt,f2ir,f3mt,f3ir)
use modmain
implicit none
! arguments
complex(8), intent(in) :: f1mt(lmmaxvr,nrmtmax,natmtot)
complex(8), intent(in) :: f1ir(ngrtot)
real(8), intent(in) :: f2mt(lmmaxvr,nrmtmax,natmtot)
real(8), intent(in) :: f2ir(ngrtot)
complex(8), intent(in) :: f3mt(lmmaxvr,nrmtmax,natmtot)
complex(8), intent(in) :: f3ir(ngrtot)
! local variables
complex(8), allocatable :: f1mt_(:,:),f3mt_(:,:)
real(8), allocatable :: f2mt_(:,:)
complex(8) zsum
integer ir,ias,lm1,lm2,lm3
complex(8) zt1,zt2
complex(8) zf1(nrmtmax)
complex(8), external :: gauntyry
allocate(f1mt_(nrmtmax,lmmaxvr))
allocate(f2mt_(nrmtmax,lmmaxvr))
allocate(f3mt_(nrmtmax,lmmaxvr))
zsum=zzero
do ir=1,ngrtot
  zsum=zsum+cfunir(ir)*dconjg(f1ir(ir))*f2ir(ir)*f3ir(ir)
end do
zsum=zsum*omega/dble(ngrtot)
! muffin-tin contribution
do ias=1,natmtot
  do lm1=1,lmmaxvr
    f1mt_(:,lm1)=f1mt(lm1,:,ias)
    f2mt_(:,lm1)=f2mt(lm1,:,ias)
    f3mt_(:,lm1)=f3mt(lm1,:,ias)
  enddo
  do lm1=1,lmmaxvr
    do lm2=1,lmmaxvr
      do lm3=1,lmmaxvr
        zt1=gauntyry(lm2l(lm1),lm2l(lm2),lm2l(lm3),&
                 lm2m(lm1),lm2m(lm2),lm2m(lm3))
        if (abs(zt1).gt.1d-12) then
          do ir=1,nrmt(ias2is(ias))
            zf1(ir)=dconjg(f1mt_(ir,lm1))*f2mt_(ir,lm2)*f3mt_(ir,lm3)*&
              spr(ir,ias2is(ias))**2
          enddo
          zt2=zzero
          do ir=1,nrmt(ias2is(ias))-1
            zt2=zt2+0.5d0*(spr(ir+1,ias2is(ias))-spr(ir,ias2is(ias)))*&
              (zf1(ir)+zf1(ir+1))
          enddo
          zsum=zsum+zt2*zt1
        endif
      enddo
    enddo
  enddo
enddo !ias
deallocate(f1mt_)
deallocate(f2mt_)
deallocate(f3mt_)
intgr_zdz=zsum
return
end function


complex(8) function lf_intgr_zdz(f1mt,f1ir,f2mt,f2ir,t,f3mt,f3ir)
use modmain
implicit none
! arguments
complex(8), intent(in) :: f1mt(lmmaxvr,nrmtmax,natmtot,ntrloc)
complex(8), intent(in) :: f1ir(ngrtot,ntrloc)
real(8), intent(in) :: f2mt(lmmaxvr,nrmtmax,natmtot,ntrloc)
real(8), intent(in) :: f2ir(ngrtot,ntrloc)
integer, intent(in) :: t(3)
complex(8), intent(in) :: f3mt(lmmaxvr,nrmtmax,natmtot,ntrloc)
complex(8), intent(in) :: f3ir(ngrtot,ntrloc)
! local variables
complex(8), allocatable :: f3mt_tmp(:,:,:)
complex(8), allocatable :: f3ir_tmp(:)
complex(8) zprod
integer it,jt,v1(3),v2(3),j,ntstep,itstep,ntloc1
integer jtloc,i,tag
logical l1
complex(8) zt1,zt2

! compute <f1|f2|f3^{T}>=\int_{-\inf}^{\inf} f1^{*}(r)f2(r)f3(r-T)dr = 
!   = \sum_{R} \int_{\Omega} f1^{*}(r+R)f2(r+R)f3(r+R-T)dr
!  f1,f3 are complex
!  f2 is real

zprod=zzero
allocate(f3mt_tmp(lmmaxvr,nrmtmax,natmtot))
allocate(f3ir_tmp(ngrtot))

j=0
ntstep=mpi_grid_map(ntr,dim_t,x=j)
do itstep=1,ntstep
  f3mt_tmp=zzero
  f3ir_tmp=zzero
  do i=0,mpi_grid_size(dim_t)-1
    ntloc1=mpi_grid_map(ntr,dim_t,x=i)
    if (itstep.le.ntloc1) then
      it=mpi_grid_map(ntr,dim_t,x=i,loc=itstep)
      v1(:)=vtl(:,it)
      v2(:)=v1(:)-t(:)
      l1=.false.
      if (v2(1).ge.tlim(1,1).and.v2(1).le.tlim(2,1).and.&
          v2(2).ge.tlim(1,2).and.v2(2).le.tlim(2,2).and.&
          v2(3).ge.tlim(1,3).and.v2(3).le.tlim(2,3)) then
        jt=ivtit(v2(1),v2(2),v2(3))
        if (jt.ne.-1) then
          l1=.true.
          jtloc=mpi_grid_map(ntr,dim_t,glob=jt,x=j)
        endif
      endif
      if (l1.and.mpi_grid_x(dim_t).eq.j.and.mpi_grid_x(dim_t).ne.i) then
        tag=(itstep*mpi_grid_size(dim_t)+i)*10
        call mpi_grid_send(f3mt(1,1,1,jtloc),lmmaxvr*nrmtmax*natmtot,&
          (/dim_t/),(/i/),tag)
        call mpi_grid_send(f3ir(1,jtloc),ngrtot,(/dim_t/),(/i/),tag+1)
      endif
      if (l1.and.mpi_grid_x(dim_t).eq.i) then
        if (j.ne.i) then
          tag=(itstep*mpi_grid_size(dim_t)+i)*10
          call mpi_grid_recieve(f3mt_tmp(1,1,1),lmmaxvr*nrmtmax*natmtot,&
            (/dim_t/),(/j/),tag)
          call mpi_grid_recieve(f3ir_tmp(1),ngrtot,(/dim_t/),(/j/),tag+1)
        else
          f3mt_tmp(:,:,:)=f3mt(:,:,:,jtloc)
          f3ir_tmp(:)=f3ir(:,jtloc)
        endif
      endif
    endif
  enddo !i
  if (itstep.le.ntrloc) then
    zt1=intgr_zdz(f1mt(1,1,1,itstep),f1ir(1,itstep),f2mt(1,1,1,itstep),&
     f2ir(1,itstep),f3mt_tmp,f3ir_tmp)
    zprod=zprod+zt1
  endif
  call mpi_grid_barrier(dims=(/dim_t/))
enddo !itstep
deallocate(f3mt_tmp,f3ir_tmp)
call mpi_grid_reduce(zprod,dims=(/dim_t/))
lf_intgr_zdz=zprod
return
end function

subroutine lf_mult_zd(alpha,f1mt,f1ir,f2mt,f2ir,f3mt,f3ir)
use modmain
implicit none
complex(8), intent(in) :: alpha
complex(8), intent(in) :: f1mt(lmmaxvr,nrmtmax,natmtot,ntrloc)
complex(8), intent(in) :: f1ir(ngrtot,ntrloc)
real(8), intent(in) :: f2mt(lmmaxvr,nrmtmax,natmtot,ntrloc)
real(8), intent(in) :: f2ir(ngrtot,ntrloc)
complex(8), intent(out) :: f3mt(lmmaxvr,nrmtmax,natmtot,ntrloc)
complex(8), intent(out) :: f3ir(ngrtot,ntrloc)
integer lm1,lm2,lm3,ir,itloc,ias
complex(8), allocatable :: f1mt_(:,:),f3mt_(:,:)
real(8), allocatable :: f2mt_(:,:)
complex(8), external :: gauntyry
complex(8) zt1

allocate(f1mt_(nrmtmax,lmmaxvr))
allocate(f2mt_(nrmtmax,lmmaxvr))
allocate(f3mt_(nrmtmax,lmmaxvr))

do itloc=1,ntrloc
  do ias=1,natmtot
    do lm1=1,lmmaxvr
      f1mt_(:,lm1)=f1mt(lm1,:,ias,itloc)
      f2mt_(:,lm1)=f2mt(lm1,:,ias,itloc)
      f3mt_(:,lm1)=zzero
    enddo
    do lm1=1,lmmaxvr
      do lm2=1,lmmaxvr
        do lm3=1,lmmaxvr
          zt1=gauntyry(lm2l(lm3),lm2l(lm2),lm2l(lm1),&
            lm2m(lm3),lm2m(lm2),lm2m(lm1))
          if (abs(zt1).gt.1d-12) then
            do ir=1,nrmt(ias2is(ias))
              f3mt_(ir,lm3)=f3mt_(ir,lm3)+f1mt_(ir,lm1)*f2mt_(ir,lm2)*zt1
            enddo
          endif
        enddo
      enddo
    enddo
    do lm3=1,lmmaxvr
      f3mt(lm3,:,ias,itloc)=alpha*f3mt_(:,lm3)
    enddo
  enddo !ias
  f3ir(:,itloc)=alpha*f1ir(:,itloc)*f2ir(:,itloc)
enddo !itloc
deallocate(f1mt_,f2mt_,f3mt_)
return
end subroutine


!complex(8) function lf_intgr_zz(f1mt,f1ir,t,f2mt,f2ir)
!use modmain
!implicit none
!complex(8), intent(in) :: f1mt(lmmaxvr,nrmtmax,natmtot,*)
!complex(8), intent(in) :: f1ir(ngrtot,*)
!integer, intent(in) :: t(3)
!complex(8), intent(in) :: f2mt(lmmaxvr,nrmtmax,natmtot,*)
!complex(8), intent(in) :: f2ir(ngrtot,*)
!complex(8), allocatable :: f2mt_tmp(:,:,:)
!complex(8), allocatable :: f2ir_tmp(:)
!
!complex(8) zprod
!integer it,jt,v1(3),v2(3),j,ntstep,ntrloc,itstep,ntloc1
!integer jtloc,i,tag
!logical l1
!complex(8), external :: zfinp_
!
!! compute <f1_0|f2_T>=\int_{-\inf}^{\inf} f1^{*}(r)f2(r-T)dr = 
!!   = \sum_{R} \int_{\Omega} f1^{*}(r+R)f2(r+R-T)dr
!
!zprod=zzero
!allocate(f2mt_tmp(lmmaxvr,nrmtmax,natmtot))
!allocate(f2ir_tmp(ngrtot))
!
!j=0
!ntstep=mpi_grid_map(ntr,dim_t,x=j)
!ntrloc=mpi_grid_map(ntr,dim_t)
!do itstep=1,ntstep
!  f2mt_tmp=zzero
!  f2ir_tmp=zzero
!  do i=0,mpi_grid_size(dim_t)-1
!    ntloc1=mpi_grid_map(ntr,dim_t,x=i)
!    if (itstep.le.ntloc1) then
!      it=mpi_grid_map(ntr,dim_t,x=i,loc=itstep)
!      v1(:)=vtl(:,it)
!      v2(:)=v1(:)-t(:)
!      l1=.false.
!      if (v2(1).ge.tlim(1,1).and.v2(1).le.tlim(2,1).and.&
!          v2(2).ge.tlim(1,2).and.v2(2).le.tlim(2,2).and.&
!          v2(3).ge.tlim(1,3).and.v2(3).le.tlim(2,3)) then
!        jt=ivtit(v2(1),v2(2),v2(3))
!        if (jt.ne.-1) then
!          l1=.true.
!          jtloc=mpi_grid_map(ntr,dim_t,glob=jt,x=j)
!        endif
!      endif
!      if (l1.and.mpi_grid_x(dim_t).eq.j.and.mpi_grid_x(dim_t).ne.i) then
!        tag=(itstep*mpi_grid_size(dim_t)+i)*10
!        call mpi_grid_send(f2mt(1,1,1,jtloc),lmmaxvr*nrmtmax*natmtot,&
!          (/dim_t/),(/i/),tag)
!        call mpi_grid_send(f2ir(1,jtloc),ngrtot,(/dim_t/),(/i/),tag+1)
!      endif
!      if (l1.and.mpi_grid_x(dim_t).eq.i) then
!        if (j.ne.i) then
!          tag=(itstep*mpi_grid_size(dim_t)+i)*10
!          call mpi_grid_recieve(f2mt_tmp(1,1,1),lmmaxvr*nrmtmax*natmtot,&
!            (/dim_t/),(/j/),tag)
!          call mpi_grid_recieve(f2ir_tmp(1),ngrtot,(/dim_t/),(/j/),tag+1)
!        else
!          f2mt_tmp(:,:,:)=f2mt(:,:,:,jtloc)
!          f2ir_tmp(:)=f2ir(:,jtloc)
!        endif
!      endif
!    endif
!  enddo !
!  if (itstep.le.ntrloc) then
!    zprod=zprod+zfinp_(.true.,f1mt(1,1,1,itstep),f2mt_tmp,f1ir(1,itstep),&
!      f2ir_tmp)
!  endif
!  call mpi_grid_barrier(dims=(/dim_t/))
!enddo
!deallocate(f2mt_tmp,f2ir_tmp)
!call mpi_grid_reduce(zprod,dims=(/dim_t/))
!lf_intgr_zz=zprod
!return
!end function






! compute f3(r)=alpha*f1(r)*f2(r)+beta*f3(r)
!subroutine lf_prod(alpha,f1mt,f1ir,f2mt,f2ir,beta,f3mt,f3ir)
!use modmain
!implicit none
!complex(8), intent(in) :: alpha
!complex(8), intent(in) :: f1mt(lmmaxvr,nrmtmax,natmtot,*)
!complex(8), intent(in) :: f1ir(ngrtot,*)
!complex(8), intent(in) :: f2mt(lmmaxvr,nrmtmax,natmtot,*)
!complex(8), intent(in) :: f2ir(ngrtot,*)
!complex(8), intent(in) :: beta
!complex(8), intent(inout) :: f3mt(lmmaxvr,nrmtmax,natmtot,*)
!complex(8), intent(inout) :: f3ir(ngrtot,*)
!integer lm1,lm2,lm3,ir,itloc,ias
!complex(8), allocatable :: f1mt_(:,:),f2mt_(:,:),f3mt_(:,:)
!real(8), external :: gaunt
!real(8) t1
!
!allocate(f1mt_(nrmtmax,lmmaxvr))
!allocate(f2mt_(nrmtmax,lmmaxvr))
!allocate(f3mt_(nrmtmax,lmmaxvr))
!
!do itloc=1,ntrloc
!  do ias=1,natmtot
!    do lm1=1,lmmaxvr
!      f1mt_(:,lm1)=f1mt(lm1,:,ias,itloc)
!      f2mt_(:,lm1)=f2mt(lm1,:,ias,itloc)
!      f3mt_(:,lm1)=beta*f3mt(lm1,:,ias,itloc)
!    enddo
!    do lm1=1,lmmaxvr
!      do lm2=1,lmmaxvr
!        do lm3=1,lmmaxvr
!          t1=gaunt(lm2l(lm3),lm2l(lm1),lm2l(lm2),&
!                   lm2m(lm3),lm2m(lm1),lm2m(lm2))
!          if (abs(t1).gt.1d-8) then
!            do ir=1,nrmt(ias2is(ias))
!              f3mt_(ir,lm3)=f3mt_(ir,lm3)+alpha*f1mt_(ir,lm1)*f2mt_(ir,lm2)*t1
!            enddo
!          endif
!        enddo
!      enddo
!    enddo
!    do lm3=1,lmmaxvr
!      f3mt(lm3,:,ias,itloc)=f3mt_(:,lm3)
!    enddo
!  enddo !ias
!  f3ir(:,itloc)=alpha*f1ir(:,itloc)*f2ir(:,itloc)+beta*f3ir(:,itloc)
!enddo !itloc
!deallocate(f1mt_,f2mt_,f3mt_)
!return
!end subroutine




!subroutine lf_sht(sht,fmt_in,fmt_out)
!use modmain
!implicit none
!character, intent(in) :: sht
!complex(8), intent(in) :: fmt_in(lmmaxvr,nrmtmax,natmtot)
!complex(8), intent(out) :: fmt_out(lmmaxvr,nrmtmax,natmtot)
!integer ias,is
!complex(8), allocatable :: f1(:,:)
!
!allocate(f1(lmmaxvr,nrmtmax))
!! backward transform (spherical harmonics to coordinates)
!if (sht.eq.'b'.or.sht.eq.'B') then
!  do ias=1,natmtot
!    is=ias2is(ias)
!    f1=fmt_in(:,:,ias)
!    call zgemm('N','N',lmmaxvr,nrmt(is),lmmaxvr,zone,zbshtvr, &
!      lmmaxvr,f1,lmmaxvr,zzero,fmt_out(:,:,ias),lmmaxvr)
!  enddo !ias
!endif
!! forward transform (spherical coordinates to harmonics)
!if (sht.eq.'f'.or.sht.eq.'F') then
!  do ias=1,natmtot
!    is=ias2is(ias)
!    f1=fmt_in(:,:,ias)
!    call zgemm('N','N',lmmaxvr,nrmt(is),lmmaxvr,zone,zfshtvr, &
!      lmmaxvr,f1,lmmaxvr,zzero,fmt_out(:,:,ias),lmmaxvr)
!  enddo !ias
!endif
!deallocate(f1)
!return
!end subroutine

subroutine lf_write(fname,fmt,fir)
use modmain
implicit none
character*(*), intent(in) :: fname
real(8), intent(in) :: fmt(lmmaxvr,nrmtmax,natmtot,ntrloc)
real(8), intent(in) :: fir(ngrtot,ntrloc)

integer it,itloc,i1,i2,i3,ir,i,is,ia,itr(3),ir0
real(8) vtrc(3),v2(3),v3(3),vrc0(3),r0
logical l1
logical, external :: vrinmt

open(160,file=trim(adjustl(fname)),status="REPLACE",form="FORMATTED")
i=0
do itloc=1,ntrloc
  it=mpi_grid_map(ntr,dim_t,loc=itloc)
  vtrc(:)=vtl(1,it)*avec(:,1)+vtl(2,it)*avec(:,2)+vtl(3,it)*avec(:,3)
  ir=0
  do i3=0,ngrid(3)-1
    v2(3)=dble(i3)/dble(ngrid(3))
    do i2=0,ngrid(2)-1
      v2(2)=dble(i2)/dble(ngrid(2))
      do i1=0,ngrid(1)-1
        v2(1)=dble(i1)/dble(ngrid(1))
        ir=ir+1
        call r3mv(avec,v2,v3)
        l1=vrinmt(v3,is,ia,itr,vrc0,ir0,r0)        
        v3(:)=v3(:)+vtrc(:)
        if (abs(fir(ir,itloc)).gt.1d-10) then
          write(160,*)v3(1),',',v3(2),',',v3(3),',',fir(ir,itloc)
          i=i+1
        endif
      enddo
    enddo
  enddo
enddo
write(*,*)'points=',i
close(160)
end subroutine


end module