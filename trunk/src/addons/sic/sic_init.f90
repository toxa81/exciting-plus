subroutine sic_init
use modmain
use mod_wannier
use mod_sic
!use mod_linresp
use mod_ws
implicit none
integer i,n,ir,j
!logical texist
!logical, external :: sic_include_cell
integer itp,lm
integer i1,i2,i3,ish,vl(3)
real(8) a,b,x0,tp(2),wsv(3,3)
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

tsicsv=.false.

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

wsv(:,1)=avec(:,1)*ngridk(1)
wsv(:,2)=avec(:,2)*ngridk(2)
wsv(:,3)=avec(:,3)*ngridk(3)
call ws_init(wsv)
sic_wan_rwsmin=1.d100
sic_wan_rwsmax=0.d0
do i=1,26
  sic_wan_rwsmin=min(sic_wan_rwsmin,sqrt(sum(ws_pts(:,i)**2/4.d0)))
enddo
do i=1,ws_nvertex
  sic_wan_rwsmax=max(sic_wan_rwsmax,sqrt(sum(ws_vertex(:,i)**2)))
enddo
if (s_rmax.gt.sic_wan_rwsmax) s_rmax=sic_wan_rwsmax
if (s_rmax.gt.sic_wan_rwsmin) then
  s_rmin=sic_wan_rwsmin
else
  s_rmin=s_rmax
endif

if (mpi_grid_root()) then
  write(*,'("[sic_init] sic_wan_rwsmin : ",G18.10)')sic_wan_rwsmin
  write(*,'("[sic_init] sic_wan_rwsmax : ",G18.10)')sic_wan_rwsmax
endif

if (allocated(sic_blochsum%vtl)) deallocate(sic_blochsum%vtl)
allocate(sic_blochsum%vtl(3,sic_maxvtl))
!sic_blochsum%ntr=0
!l2=.true.
!ish=0
!do while (l2)
!  l1=.false.
!  do i1=-ish,ish
!    do i2=-ish,ish
!      do i3=-ish,ish
!        if (abs(i1).eq.ish.or.abs(i2).eq.ish.or.abs(i3).eq.ish) then
!          vl=(/i1,i2,i3/)
!          if (sic_include_cell(vl)) then
!            l1=.true.
!            sic_blochsum%ntr=sic_blochsum%ntr+1
!            if (sic_blochsum%ntr.gt.sic_maxvtl) then
!              write(*,'("Error(sic_init) : sic_maxvtl is too small")')
!              call pstop
!            endif
!            sic_blochsum%vtl(:,sic_blochsum%ntr)=vl
!          endif
!        endif
!      enddo
!    enddo
!  enddo !i1
!  if (l1) then
!    ish=ish+1
!  else
!    l2=.false.
!  endif
!enddo

sic_blochsum%ntr=0
do i1=0,ngridk(1)-1
  do i2=0,ngridk(2)-1
    do i3=0,ngridk(3)-1
      sic_blochsum%ntr=sic_blochsum%ntr+1
      sic_blochsum%vtl(:,sic_blochsum%ntr)=(/i1,i2,i3/)
    enddo
  enddo
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

if (allocated(s_wlm)) deallocate(s_wlm)
allocate(s_wlm(lmmaxwan,s_nr,nspinor,sic_wantran%nwan))
s_wlm=zzero
if (allocated(s_wvlm)) deallocate(s_wvlm)
allocate(s_wvlm(lmmaxwan,s_nr,nspinor,sic_wantran%nwan))
s_wvlm=zzero
if (allocated(sic_vme)) deallocate(sic_vme)
allocate(sic_vme(sic_wantran%nwt))
sic_vme=zzero
if (allocated(sic_wb)) deallocate(sic_wb)
if (allocated(sic_wvb)) deallocate(sic_wvb)
if (tsveqn) then
  allocate(sic_wb(sic_wantran%nwan,nstfv,nspinor,nkptloc))
  allocate(sic_wvb(sic_wantran%nwan,nstfv,nspinor,nkptloc))
else
  if (tsicsv) then
    allocate(sic_wb(sic_wantran%nwan,nstsv,1,nkptloc))
    allocate(sic_wvb(sic_wantran%nwan,nstsv,1,nkptloc))  
  else
    allocate(sic_wb(sic_wantran%nwan,nmatmax,nspinor,nkptloc))
    allocate(sic_wvb(sic_wantran%nwan,nmatmax,nspinor,nkptloc))
  endif
endif
sic_wb=zzero
sic_wvb=zzero
if (.not.allocated(sic_wan_e0)) then
  allocate(sic_wan_e0(nwantot))
  sic_wan_e0=0.d0
endif
if (.not.allocated(sic_wan_umtrx)) then
  allocate(sic_wan_umtrx(nwantot,nwantot,nkptnrloc))
  sic_wan_umtrx=zzero
  do i=1,nwantot
    sic_wan_umtrx(i,i,:)=zone
  enddo
endif
if (allocated(sic_wan_h0k)) deallocate(sic_wan_h0k)
allocate(sic_wan_h0k(sic_wantran%nwan,sic_wantran%nwan,nkptloc))
sic_wan_h0k=zzero

if (allocated(s_wkmt)) deallocate(s_wkmt)
allocate(s_wkmt(mt_ntp,nrmtmax,natmtot,nspinor,sic_wantran%nwan,nkptloc))
if (allocated(s_wkit)) deallocate(s_wkit)
allocate(s_wkit(ngkmax,nspinor,sic_wantran%nwan,nkptloc))
if (allocated(s_wvkmt)) deallocate(s_wvkmt)
allocate(s_wvkmt(mt_ntp,nrmtmax,natmtot,nspinor,sic_wantran%nwan,nkptloc))
if (allocated(s_wvkit)) deallocate(s_wvkit)
allocate(s_wvkit(ngkmax,nspinor,sic_wantran%nwan,nkptloc))
if (.not.tsveqn) then
  if (allocated(s_wkmtlm)) deallocate(s_wkmtlm)
  allocate(s_wkmtlm(nrmtmax,lmmaxapw,natmtot,nspinor,sic_wantran%nwan,nkptloc))
  if (allocated(s_wvkmtlm)) deallocate(s_wvkmtlm)
  allocate(s_wvkmtlm(nrmtmax,lmmaxapw,natmtot,nspinor,sic_wantran%nwan,nkptloc))
endif
if (mpi_grid_root()) then
! size of sperical arrays
  a=2*16.d0*lmmaxwan*s_nr*nspinor*sic_wantran%nwan/1024/1024
! size of Bloch sums
  a=a+2*16.d0*mt_ntp*nrmtmax*natmtot*nspinor*sic_wantran%nwan*nkptloc/1024/1024
  a=a+2*16.d0*ngkmax*nspinor*sic_wantran%nwan*nkptloc/1024/1204
  write(*,'("[sic_init] Memory usage (Mb) : ",F12.4)')a
endif
return
end subroutine

!logical function sic_include_cell(vl)
!use modmain
!use mod_sic
!implicit none
!integer, intent(in) :: vl(3)
!logical l1
!integer n,jas,ias,ir,j
!real(8) vt(3),v1(3)
!vt(:)=vl(1)*avec(:,1)+vl(2)*avec(:,2)+vl(3)*avec(:,3)
!l1=.false.
!do j=1,sic_wantran%nwan
!  n=sic_wantran%iwan(j)
!  do ias=1,natmtot
!    v1(:)=atposc(:,ias2ia(ias),ias2is(ias))+vt(:)- wanpos(:,n)
!    if (sqrt(sum(v1(:)**2)).le.(sic_wan_cutoff+rmt(ias2is(ias)))) l1=.true.
!  enddo
!enddo
!sic_include_cell=l1
!return
!end
!

