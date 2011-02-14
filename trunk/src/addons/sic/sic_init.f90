subroutine sic_init
use modmain
use mod_wannier
use mod_sic
use mod_linresp
use mod_wannier
implicit none
integer i,n,ir,j
logical texist
logical, external :: vrinmt,sic_include_cell
integer itp,lm
integer i1,i2,i3,ish,vl(3)
real(8) a,b,x0,tp(2)
logical l1,l2
!
tevecsv=.true.
lmmaxwan=(lmaxwan+1)**2

if (allocated(sic_apply)) deallocate(sic_apply)
allocate(sic_apply(nwantot))
if (allocated(sicw)) then
  sic_apply=0
  do n=1,nwantot
    i=wan_info(6,n)
    j=wan_info(7,n)
    sic_apply(n)=sicw(j,i)
  enddo
else
  sic_apply=2
endif
call deletewantran(sic_wantran)
! get Wannier transitions
call genwantran(sic_wantran,-0.d0,sic_me_cutoff,allwt=.true.,waninc=sic_apply)

if (allocated(sic_orbitals%vtl)) deallocate(sic_orbitals%vtl)
allocate(sic_orbitals%vtl(3,sic_maxvtl))
sic_orbitals%vtl=-1000000
sic_orbitals%ntr=0
l2=.true.
ish=0
do while (l2)
  l1=.false.
  do i1=-ish,ish
    do i2=-ish,ish
      do i3=-ish,ish
        if (abs(i1).eq.ish.or.abs(i2).eq.ish.or.abs(i3).eq.ish) then
          vl=(/i1,i2,i3/)
          if (sic_include_cell(vl)) then
            l1=.true.            
            sic_orbitals%ntr=sic_orbitals%ntr+1
            if (sic_orbitals%ntr.gt.sic_maxvtl) then
              write(*,'("Error(sic_init) : sic_maxvtl is too small")')
              call pstop
            endif
            sic_orbitals%vtl(:,sic_orbitals%ntr)=vl
          endif
        endif
      enddo
    enddo
  enddo !i1
  if (l1) then
    ish=ish+1
  else
    l2=.false.
  endif
enddo
! Cartesian coordinates of translation vectors
if (allocated(sic_orbitals%vtc)) deallocate(sic_orbitals%vtc)
allocate(sic_orbitals%vtc(3,sic_orbitals%ntr))
do i=1,sic_orbitals%ntr
  sic_orbitals%vtc(:,i)=sic_orbitals%vtl(1,i)*avec(:,1)+&
                        sic_orbitals%vtl(2,i)*avec(:,2)+&
                        sic_orbitals%vtl(3,i)*avec(:,3)
end do
!! find translation limits
!sic_orbitals%tlim=0
!do i=1,3
!  sic_orbitals%tlim(1,i)=minval(sic_orbitals%vtl(i,1:sic_orbitals%ntr))
!  sic_orbitals%tlim(2,i)=maxval(sic_orbitals%vtl(i,1:sic_orbitals%ntr))
!enddo
!! find mapping from translations to linear index
!if (allocated(sic_orbitals%ivtit)) deallocate(sic_orbitals%ivtit)
!allocate(sic_orbitals%ivtit(sic_orbitals%tlim(1,1):sic_orbitals%tlim(2,1),&
!                            sic_orbitals%tlim(1,2):sic_orbitals%tlim(2,2),&
!                            sic_orbitals%tlim(1,3):sic_orbitals%tlim(2,3)))
!sic_orbitals%ivtit=-1
!do i=1,sic_orbitals%ntr
!  sic_orbitals%ivtit(sic_orbitals%vtl(1,i),sic_orbitals%vtl(2,i),&
!    sic_orbitals%vtl(3,i))=i
!enddo

call sic_genrmesh

if (allocated(s_tpw)) deallocate(s_tpw)
allocate(s_tpw(s_ntp))
if (allocated(s_rlmf)) deallocate(s_rlmf)
allocate(s_rlmf(lmmaxwan,s_ntp))
if (allocated(s_ylmf)) deallocate(s_ylmf)
allocate(s_ylmf(lmmaxwan,s_ntp))
if (allocated(s_rlmb)) deallocate(s_rlmb)
allocate(s_rlmb(s_ntp,lmmaxwan))
if (allocated(s_ylmb)) deallocate(s_ylmb)
allocate(s_ylmb(s_ntp,lmmaxwan))
if (allocated(s_spx)) deallocate(s_spx)
allocate(s_spx(3,s_ntp))
if (allocated(s_wanlm)) deallocate(s_wanlm)
allocate(s_wanlm(lmmaxwan,s_nr,nspinor,sic_wantran%nwan))
if (allocated(s_wvlm)) deallocate(s_wvlm)
allocate(s_wvlm(lmmaxwan,s_nr,nspinor,sic_wantran%nwan))
if (allocated(vwanme)) deallocate(vwanme)
allocate(vwanme(sic_wantran%nwt))
vwanme=zzero
if (allocated(sic_wb)) deallocate(sic_wb)
allocate(sic_wb(sic_wantran%nwan,nstfv,nspinor,nkptloc))
sic_wb=zzero
if (allocated(sic_wvb)) deallocate(sic_wvb)
allocate(sic_wvb(sic_wantran%nwan,nstfv,nspinor,nkptloc))
sic_wvb=zzero
if (allocated(sic_wann_e0)) deallocate(sic_wann_e0)
allocate(sic_wann_e0(nwantot))
sic_wann_e0=0.d0
if (allocated(sic_wann_h0k)) deallocate(sic_wann_h0k)
allocate(sic_wann_h0k(sic_wantran%nwan,sic_wantran%nwan,nkptloc))
sic_wann_h0k=zzero
if (allocated(sic_wgk)) deallocate(sic_wgk)
allocate(sic_wgk(ngkmax,sic_wantran%nwan,nspinor,nkptloc))
sic_wgk=zzero
if (allocated(sic_wvgk)) deallocate(sic_wvgk)
allocate(sic_wvgk(ngkmax,sic_wantran%nwan,nspinor,nkptloc))
sic_wgk=zzero
if (allocated(sic_wuy)) deallocate(sic_wuy)
allocate(sic_wuy(lmmaxvr,nufrmax,natmtot,sic_wantran%nwan,nspinor,nkptloc))
sic_wuy=zzero
if (allocated(sic_wvuy)) deallocate(sic_wvuy)
allocate(sic_wvuy(lmmaxvr,nufrmax,natmtot,sic_wantran%nwan,nspinor,nkptloc))
sic_wvuy=zzero

if (allocated(s_wankmt)) deallocate(s_wankmt)
allocate(s_wankmt(lmmaxvr,nrmtmax,natmtot,nspinor,sic_wantran%nwan,nkptloc))
if (allocated(s_wankir)) deallocate(s_wankir)
allocate(s_wankir(ngrtot,nspinor,sic_wantran%nwan,nkptloc))
if (allocated(s_wvkmt)) deallocate(s_wvkmt)
allocate(s_wvkmt(lmmaxvr,nrmtmax,natmtot,nspinor,sic_wantran%nwan,nkptloc))
if (allocated(s_wvkir)) deallocate(s_wvkir)
allocate(s_wvkir(ngrtot,nspinor,sic_wantran%nwan,nkptloc))

inquire(file="SIC_WANN_E0.OUT",exist=texist)
if (texist) then
  open(170,file="SIC_WANN_E0.OUT",form="FORMATTED",status="OLD")
  do n=1,nwantot
    read(170,*,err=20)sic_wann_e0(n)
  enddo
  close(170)
  goto 30
20 sic_wann_e0=0.d0
endif
30 continue
! Lebedev-Laikov mesh
call leblaik(s_ntp,s_spx,s_tpw)
! get (theta,phi) of each spx vector and generate spherical harmonics
do itp=1,s_ntp
  s_tpw(itp)=s_tpw(itp)*fourpi
  call sphcrd(s_spx(:,itp),a,tp)
  call genrlm(lmaxwan,tp,s_rlmf(1,itp))
  call genylm(lmaxwan,tp,s_ylmf(1,itp))
  do lm=1,lmmaxwan
    s_rlmb(itp,lm)=s_rlmf(lm,itp)*s_tpw(itp)
    s_ylmb(itp,lm)=dconjg(s_ylmf(lm,itp))*s_tpw(itp) 
  enddo
enddo
!do i=4,10
!  call test1(6,i)
!  call test1(14,i)
!  call test1(26,i)
!  call test1(38,i)
!  call test1(50,i)
!  call test1(74,i)
!  call test1(110,i)
!  call test1(146,i) 
!  call test1(170,i)
!  call test1(194,i)
!  call test1(230,i)
!  call test1(266,i)
!  call test1(302,i)
!  call test1(350,i)
!  call test1(434,i)
!  call test1(590,i)
!  call test1(770,i)
!  call test1(974,i)
!  call test1(1202,i)
!enddo
!call bstop

if (mpi_grid_root()) then
! size of sperical arrays
  a=2*16.d0*lmmaxwan*s_nr*nspinor*sic_wantran%nwan/1024/1024
!! size of Bloch sums
!  a=a+2*16.d0*lmmaxvr*nrmtmax*natmtot*nspinor*sic_wantran%nwan*nkptloc/1024/1024
!  a=a+2*16.d0*ngrtot*nspinor*sic_wantran%nwan*nkptloc/1024/1204
  write(*,'("[sic_init] Memory usage (Mb) : ",F12.4)')a
endif
return
end

logical function sic_include_cell(vl)
use modmain
use mod_sic
implicit none
integer, intent(in) :: vl(3)
logical l1
integer n,ias,jas,ir,j
real(8) vt(3),v1(3)
vt(:)=vl(1)*avec(:,1)+vl(2)*avec(:,2)+vl(3)*avec(:,3)
l1=.false.
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  ias=wan_info(1,n)
  do jas=1,natmtot
    v1(:)=atposc(:,ias2ia(jas),ias2is(jas))+vt(:)-&
      atposc(:,ias2ia(ias),ias2is(ias))
    if (sqrt(sum(v1(:)**2)).le.sic_wan_cutoff) l1=.true.
  enddo
  do ir=1,ngrtot
    v1(:)=vgrc(:,ir)+vt(:)-atposc(:,ias2ia(ias),ias2is(ias))
    if (sqrt(sum(v1(:)**2)).le.sic_wan_cutoff) l1=.true.
  enddo
enddo
sic_include_cell=l1
return
end

!
!subroutine test1(ntp,lmax)
!use modmain
!implicit none
!integer, intent(in) :: ntp
!integer, intent(in) :: lmax
!integer itp,lm,lmmax,itp1
!real(8) a,tp(2)
!complex(8) z1
!real(8), allocatable :: spx(:,:)
!real(8), allocatable :: tpw(:)
!complex(8), allocatable :: ylmf(:,:)
!complex(8), allocatable :: ylmb(:,:)
!
!lmmax=(lmax+1)**2
!
!allocate(spx(3,ntp))
!allocate(tpw(ntp))
!allocate(ylmf(lmmax,ntp))
!allocate(ylmb(ntp,lmmax))
!! Lebedev-Laikov mesh
!call leblaik(ntp,spx,tpw)
!! get (theta,phi) of each spx vector and generate spherical harmonics
!do itp=1,ntp                   
!  tpw(itp)=tpw(itp)*fourpi
!  call sphcrd(spx(1,itp),a,tp)
!  call genylm(lmax,tp,ylmf(1,itp))
!  do lm=1,lmmax
!    ylmb(itp,lm)=dconjg(ylmf(lm,itp))*tpw(itp) 
!  enddo  
!enddo 
!a=0.d0
!do itp=1,ntp
!  do itp1=1,ntp
!    z1=zzero
!    do lm=1,lmmax
!      z1=z1+ylmb(itp,lm)*ylmf(lm,itp1)
!    enddo
!    !if (itp.eq.itp1) z1=z1-zone
!    !a=max(a,abs(z1))
!    if (itp.ne.itp1) a=max(a,abs(z1))
!  enddo
!enddo
!if (mpi_grid_root()) then
!  write(*,'("[test1] l = ",I4," ntp = ",I4)')lmax,ntp
!  write(*,'("[test1] Lebedev quadrature completeness error : ",G18.10)')a
!  write(*,*)
!endif
!deallocate(spx,tpw,ylmf,ylmb)
!return
!end

