subroutine sic_init
use modmain
use mod_wannier
use mod_sic
use mod_ws
implicit none
integer i,n,j,i1,i2,i3,ivl(3)
real(8) a,wsv(3,3),vl(3),vc(3)
!
if (.not.wannier) then
  write(*,*)
  write(*,'("Error(sic_init): no Wannier functions")')
  write(*,*)
  call pstop
endif
tevecsv=.true.
lmmaxwan=(lmaxwan+1)**2

sic_bottom_energy=-20.3d0

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
  sic_wan_rwsmin=min(sic_wan_rwsmin,sqrt(sum(ws_nnpts(:,i)**2/4.d0)))
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
call sic_genrmesh
call sic_gensmesh
! translation vectors for Bloch-sums
sic_ntr=ngridk(1)*ngridk(2)*ngridk(3)
if (allocated(sic_vtl)) deallocate(sic_vtl)
allocate(sic_vtl(3,sic_ntr))
if (allocated(sic_vtc)) deallocate(sic_vtc)
allocate(sic_vtc(3,sic_ntr))
i=0
do i1=0,ngridk(1)-1
  do i2=0,ngridk(2)-1
    do i3=0,ngridk(3)-1
      ivl(:)=(/i1,i2,i3/)
      vl(:)=dble(ivl(:))
      call r3mv(avec,vl,vc)
      i=i+1
      sic_vtl(:,i)=ivl(:)
      sic_vtc(:,i)=vc(:)
    enddo
  enddo
enddo
! main arrays
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
  allocate(sic_wb(sic_wantran%nwan,nmatmax,nspinor,nkptloc))
  allocate(sic_wvb(sic_wantran%nwan,nmatmax,nspinor,nkptloc))
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
allocate(s_wkmt(nrmtmax,lmmaxapw,natmtot,nspinor,sic_wantran%nwan,nkptloc))
if (allocated(s_wkit)) deallocate(s_wkit)
allocate(s_wkit(ngkmax,nspinor,sic_wantran%nwan,nkptloc))
if (allocated(s_wvkmt)) deallocate(s_wvkmt)
allocate(s_wvkmt(nrmtmax,lmmaxapw,natmtot,nspinor,sic_wantran%nwan,nkptloc))
if (allocated(s_wvkit)) deallocate(s_wvkit)
allocate(s_wvkit(ngkmax,nspinor,sic_wantran%nwan,nkptloc))
if (mpi_grid_root()) then
! size of spherical arrays
  a=2*16.d0*lmmaxwan*s_nr*nspinor*sic_wantran%nwan/1024/1024
! size of Bloch sums
  a=a+2*16.d0*lmmaxapw*nrmtmax*natmtot*nspinor*sic_wantran%nwan*nkptloc/1024/1024
  a=a+2*16.d0*ngkmax*nspinor*sic_wantran%nwan*nkptloc/1024/1204
  write(*,'("[sic_init] Memory usage (Mb) : ",F12.4)')a
endif
return
end subroutine
