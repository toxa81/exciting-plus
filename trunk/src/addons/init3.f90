subroutine init3
use modmain
use mod_wannier
implicit none
integer ia,is,l,m,ir,i1,i2,i3,i
integer l1,l2,l3,m1,m2,m3,lm1,lm2,lm3
real(8) vl(3)
complex(8), external :: gauntyry
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
if (allocated(sv_gntyry)) deallocate(sv_gntyry)
allocate(sv_gntyry(lmmaxvr,lmmaxvr,lmmaxvr))
do l1=0,lmaxvr
  do m1=-l1,l1
    lm1=idxlm(l1,m1)
    do l2=0,lmaxvr
      do m2=-l2,l2
        lm2=idxlm(l2,m2)
        do l3=0,lmaxvr
          do m3=-l3,l3
            lm3=idxlm(l3,m3)
            sv_gntyry(lm3,lm2,lm1)=gauntyry(l2,l3,l1,m2,m3,m1)
          enddo
        enddo
      enddo
    enddo
  enddo
enddo
if (allocated(sv_ubu)) deallocate(sv_ubu)
allocate(sv_ubu(lmmaxvr,nlufrmax,nlufrmax,natmtot,ndmag))
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
