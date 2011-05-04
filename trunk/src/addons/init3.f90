subroutine init3
use modmain
use mod_wannier
implicit none
integer ia,is,l,m,ir,i1,i2,i3,i,n
logical l1(maxspst,0:lmaxapw)
real(8) vl(3)
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
call getnufr(lmaxapw)
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
!if (allocated(spnl)) deallocate(spnl)
!allocate(spnl(0:maxlapw,maxspecies))
!spnl=0
!do l=0,lmaxapw
!  spnl(l,:)=l+1
!enddo
!do is=1,nspecies
!  l1=.false.
!  do i=1,spnst(is)
!    n=spn(i,is)
!    l=spl(i,is)
!    if (spcore(i,is).and..not.l1(n,l)) then
!      spnl(l,is)=spnl(l,is)+1
!      l1(n,l)=.true.
!    endif
!  enddo
!enddo
if (wannier) call wann_init
if (sic) call sic_init
if (debug_level.ge.4) then
  fdbgout=999
  write(fdbgname,'("iproc_",I7.7,"__debug.txt")')iproc
  open(fdbgout,file=trim(adjustl(fdbgname)),form="FORMATTED",status="REPLACE")
  do i=1,mpi_grid_nd
    write(fdbgout,'("x(",I2,") : ",I8)')i,mpi_grid_dim_pos(i)
  enddo
  write(fdbgout,'("task : ",I4)')task
  close(fdbgout)
endif
!if (mpi_grid_root()) call srclog
!if (mpi_grid_root()) call print_info 
return
end
