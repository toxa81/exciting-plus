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
if (.not.wannier) then
  write(*,*)
  write(*,'("Error(sic_init): no Wannier functions")')
  write(*,*)
  call pstop
endif
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

if (allocated(sic_blochsum%vtl)) deallocate(sic_blochsum%vtl)
allocate(sic_blochsum%vtl(3,sic_maxvtl))
sic_blochsum%ntr=0
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
            sic_blochsum%ntr=sic_blochsum%ntr+1
            if (sic_blochsum%ntr.gt.sic_maxvtl) then
              write(*,'("Error(sic_init) : sic_maxvtl is too small")')
              call pstop
            endif
            sic_blochsum%vtl(:,sic_blochsum%ntr)=vl
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
if (allocated(sic_blochsum%vtc)) deallocate(sic_blochsum%vtc)
allocate(sic_blochsum%vtc(3,sic_blochsum%ntr))
do i=1,sic_blochsum%ntr
  sic_blochsum%vtc(:,i)=sic_blochsum%vtl(1,i)*avec(:,1)+&
                        sic_blochsum%vtl(2,i)*avec(:,2)+&
                        sic_blochsum%vtl(3,i)*avec(:,3)
end do

call sic_genrmesh
call sic_gensmesh

if (allocated(s_wanlm)) deallocate(s_wanlm)
allocate(s_wanlm(lmmaxwan,s_nr,nspinor,sic_wantran%nwan))
s_wanlm=zzero
if (allocated(s_wvlm)) deallocate(s_wvlm)
allocate(s_wvlm(lmmaxwan,s_nr,nspinor,sic_wantran%nwan))
s_wvlm=zzero
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

if (allocated(s_wankmt)) deallocate(s_wankmt)
allocate(s_wankmt(mt_ntp,nrmtmax,natmtot,nspinor,sic_wantran%nwan,nkptloc))
if (allocated(s_wankir)) deallocate(s_wankir)
allocate(s_wankir(ngkmax,nspinor,sic_wantran%nwan,nkptloc))
if (allocated(s_wvkmt)) deallocate(s_wvkmt)
allocate(s_wvkmt(mt_ntp,nrmtmax,natmtot,nspinor,sic_wantran%nwan,nkptloc))
if (allocated(s_wvkir)) deallocate(s_wvkir)
allocate(s_wvkir(ngkmax,nspinor,sic_wantran%nwan,nkptloc))

!inquire(file="SIC_WANN_E0.OUT",exist=texist)
!if (texist) then
!  open(170,file="SIC_WANN_E0.OUT",form="FORMATTED",status="OLD")
!  do n=1,nwantot
!    read(170,*,err=20,end=20)sic_wann_e0(n)
!  enddo
!  close(170)
!  goto 30
!20 sic_wann_e0=0.d0
!  close(170)
!endif
!30 continue
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
integer n,jas,ias,ir,j
real(8) vt(3),v1(3)
vt(:)=vl(1)*avec(:,1)+vl(2)*avec(:,2)+vl(3)*avec(:,3)
l1=.false.
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  do ias=1,natmtot
    v1(:)=atposc(:,ias2ia(ias),ias2is(ias))+vt(:)- wanpos(:,n)
    if (sqrt(sum(v1(:)**2)).le.(sic_wan_cutoff+rmt(ias2is(ias)))) l1=.true.
  enddo
enddo
sic_include_cell=l1
return
end


