subroutine init3
use modmain
use mod_wannier
use mod_sic
use mod_rs
use mod_libapw
implicit none
integer ia,is,l,m,ir,i1,i2,i3,i
integer itp
real(8) vl(3),x1,x2,x3,tp(2),a
!
if (allocated(rylm)) deallocate(rylm)
allocate(rylm(16,16))
if (allocated(yrlm)) deallocate(yrlm)
allocate(yrlm(16,16))
if (allocated(rylm_lps)) deallocate(rylm_lps)
allocate(rylm_lps(16,16,natmtot))
if (allocated(yrlm_lps)) deallocate(yrlm_lps)
allocate(yrlm_lps(16,16,natmtot))
call genshmat
if (allocated(ias2is)) deallocate(ias2is)
allocate(ias2is(natmtot))
if (allocated(ias2ia)) deallocate(ias2ia)
allocate(ias2ia(natmtot))
do is=1,nspecies
  do ia=1,natoms(is)
    ias2is(idxas(ia,is))=is
    ias2ia(idxas(ia,is))=ia
  end do
end do
call getatmcls
if (allocated(nufr)) deallocate(nufr)
allocate(nufr(0:lmaxapw,nspecies))
if (allocated(nlufr)) deallocate(nlufr)
allocate(nlufr(nspecies))
call getnufr
if (allocated(ufr)) deallocate(ufr)
allocate(ufr(nrmtmax,0:lmaxapw,nufrmax,natmcls))
if (allocated(ufrp)) deallocate(ufrp)
allocate(ufrp(0:lmaxapw,nufrmax,nufrmax,natmcls))
if (allocated(lm2l)) deallocate(lm2l)
allocate(lm2l(51**2))
if (allocated(lm2m)) deallocate(lm2m)
allocate(lm2m(51**2))
do l=0,50
  do m=-l,l
    lm2l(idxlm(l,m))=l
    lm2m(idxlm(l,m))=m
  end do
end do
if (allocated(vgrc)) deallocate(vgrc)
allocate(vgrc(3,ngrtot))
! Cartesian coordinates of FFT grid r-vectors
ir=0
do i3=0,ngrid(3)-1
  vl(3)=dble(i3)/dble(ngrid(3))
  do i2=0,ngrid(2)-1
    vl(2)=dble(i2)/dble(ngrid(2))
    do i1=0,ngrid(1)-1
      vl(1)=dble(i1)/dble(ngrid(1))
      ir=ir+1
      call r3mv(avec,vl,vgrc(1,ir))
    enddo
  enddo
enddo
! generate Lebedev mesh for muffin-tins
if (allocated(mt_spx)) deallocate(mt_spx)
allocate(mt_spx(3,mt_ntp))
if (allocated(mt_tpw)) deallocate(mt_tpw)
allocate(mt_tpw(mt_ntp))
if (allocated(mt_ylmf)) deallocate(mt_ylmf)
allocate(mt_ylmf(lmmaxapw,mt_ntp))
if (allocated(mt_ylmb)) deallocate(mt_ylmb)
allocate(mt_ylmb(mt_ntp,lmmaxapw))
call leblaik(mt_ntp,mt_spx,mt_tpw)
do itp=1,mt_ntp
  mt_tpw(itp)=mt_tpw(itp)*fourpi
  call sphcrd(mt_spx(:,itp),a,tp)
  call genylm(lmaxapw,tp,mt_ylmf(1,itp)) 
  mt_ylmb(itp,:)=dconjg(mt_ylmf(:,itp))*mt_tpw(itp)
enddo
! radial weights for muffin-tins
if (allocated(mt_rw)) deallocate(mt_rw)
allocate(mt_rw(nrmtmax,nspecies))
do is=1,nspecies
  x2=spr(1,is)
  x3=spr(2,is)
  mt_rw(1,is)=-((x2-x3)*(3*x2**2+2*x2*x3+x3**2))/12.d0
  do ir=2,nrmt(is)-1
    x1=spr(ir-1,is)
    x2=spr(ir,is)
    x3=spr(ir+1,is)
    mt_rw(ir,is)=-((x1-x3)*(x1**2+x2**2+x2*x3+x3**2+x1*(x2+x3)))/12.d0
  enddo
  x1=spr(nrmt(is)-1,is)
  x2=spr(nrmt(is),is)
  mt_rw(nrmt(is),is)=-((x1-x2)*(x1**2+2*x1*x2+3*x2**2))/12.d0
  if (abs(sum(mt_rw(1:nrmt(is),is))-(rmt(is)**3)/3).gt.1d-10) then
    write(*,'("Error(init3): wrong weight for is : ",I4)')is
    call pstop
  endif
enddo 
if (texactrho) then
  if (allocated(sv_ubu)) deallocate(sv_ubu)
  allocate(sv_ubu(lmmaxvr,nlufrmax,nlufrmax,natmtot,ndmag))
  if (allocated(rhomagmt)) deallocate(rhomagmt)
  allocate(rhomagmt(nlufrmax,nlufrmax,lmmaxvr,natmtot,nspinor))
  if (allocated(rhomagit)) deallocate(rhomagit)
  allocate(rhomagit(ngrtot,nspinor))
endif
if (.not.tsveqn.and.spinpol) then
  if (allocated(baa)) deallocate(baa)
  allocate(baa(lmmaxvr,apwordmax,0:lmaxapw,apwordmax,0:lmaxapw,natmtot,ndmag))
  if (allocated(bloa)) deallocate(bloa)
  allocate(bloa(lmmaxvr,nlomax,apwordmax,0:lmaxapw,natmtot,ndmag))
  if (allocated(blolo)) deallocate(blolo)
  allocate(blolo(lmmaxvr,nlomax,nlomax,natmtot,ndmag))
  if (allocated(beffig)) deallocate(beffig)
  allocate(beffig(ngvec,ndmag))
endif
if (wannier) call wann_init
if (sic) call sic_init
call rs_init
! init xml_info variables
if (.not.allocated(xml_info%magmom)) allocate(xml_info%magmom(natmtot))
if (wannier) then
  if (.not.allocated(xml_info%wan_spread)) allocate(xml_info%wan_spread(nwantot))
endif
! init debug output file
if (debug_level.ge.4.and..not.tdbgout_init) then
  fdbgout=999
  write(fdbgname,'("iproc_",I7.7,"__debug.txt")')iproc
  open(fdbgout,file=trim(adjustl(fdbgname)),form="FORMATTED",&
    status="REPLACE")
  close(fdbgout)
  tdbgout_init=.true.
endif
if (debug_level.ge.4) then
  call dbg_open_file
  do i=1,mpi_grid_nd
    write(fdbgout,'("x(",I2,") : ",I8)')i,mpi_grid_dim_pos(i)
  enddo
  write(fdbgout,'("task : ",I4)')task
  call dbg_close_file
endif
!if (mpi_grid_root()) call srclog
!if (mpi_grid_root()) call print_info
call libapw_init
return
end
