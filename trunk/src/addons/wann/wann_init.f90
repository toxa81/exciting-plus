subroutine wann_init
use modmain
use modldapu
use mod_wannier
implicit none

integer i,j,n,lm,ispn,iwgrp,itype,iatom
integer, allocatable :: iwann_tmp(:,:)

if (nspnfv.eq.2) then
  write(*,*)
  write(*,'("Error(wann_init): can''t do Wannier + spin spirals")')
  write(*,*)
  call pstop
endif

nwantot=0
do i=1,wann_natom
  iwgrp=wann_iprj(2,i)
  nwantot=nwantot+wann_norb(iwgrp)
enddo
if (allocated(wan_info)) deallocate(wan_info)
allocate(wan_info(7,nwantot))
n=0
do i=1,wann_natom
  iatom=wann_iprj(1,i)
  iwgrp=wann_iprj(2,i)
  do j=1,wann_norb(iwgrp)
    n=n+1
    lm=wann_iorb(1,j,iwgrp)
    if (lm.lt.1.or.lm.gt.16) then
      write(*,*)
      write(*,'("Error(wann_init) : lm- value is out of range ")')
      write(*,*)
      call pstop
    endif
    ispn=wann_iorb(2,j,iwgrp)
    if ((spinpol.and.(ispn.lt.1.or.ispn.gt.2)).or.(.not.spinpol.and.ispn.ne.1)) then
      write(*,*)
      write(*,'("Error(wann_init) : spin projection value is out of range ")')
      write(*,*)
      call pstop
    endif
    itype=wann_iorb(3,j,iwgrp)
    wan_info(1,n)=iatom
    wan_info(2,n)=lm
    wan_info(3,n)=ispn
    wan_info(4,n)=itype
    wan_info(5,n)=iwgrp
    wan_info(6,n)=i
    wan_info(7,n)=j
  enddo
enddo !i
! sort WF by spin (for collinear case)
!if (.not.ncmag) then
!  allocate(iwann_tmp(7,nwantot))
!  n=0
!  do j=1,nwantot
!    if (wan_info(3,j).eq.1) then
!      n=n+1
!      iwann_tmp(:,n)=wan_info(:,j)
!    endif
!  enddo
!  do j=1,nwantot
!    if (wan_info(3,j).eq.2) then
!      n=n+1
!      iwann_tmp(:,n)=wan_info(:,j)
!    endif
!  enddo
!  wan_info=iwann_tmp
!  deallocate(iwann_tmp)
!endif

if (allocated(nwannias)) deallocate(nwannias)
allocate(nwannias(natmtot))
nwannias=0
if (mpi_grid_root()) then     
  open(100,file='WANNIER.OUT',form='formatted',status='replace')
  write(100,'("number of WF atoms : ", I4)')wann_natom
  write(100,'("number of WF orbital groups : ", I4)')wann_norbgrp
  write(100,'("number of WF types : ", I4)')wann_ntype
  write(100,*)
  write(100,'("total number of WF : ", I4)')nwantot
  write(100,*)
endif
do n=1,nwantot
  iatom=wan_info(1,n)
  lm=wan_info(2,n)
  ispn=wan_info(3,n)
  itype=wan_info(4,n)
  nwannias(iatom)=nwannias(iatom)+1
  if (mpi_grid_root()) then     
    write(100,'("  wf : ",I4)')n
    write(100,'("    type : ",I4)')itype
    write(100,'("    pure spinor orbital for projection : ")')
    write(100,'("      atom : ",I4)')iatom
    write(100,'("      l,m  : ",2I4)')lm2l(lm),lm-lm2l(lm)**2
    write(100,'("      ispn : ",I4)')ispn
    write(100,'("  interval : [",F8.4,",",F8.4,"]")')wann_eint(:,itype)
    write(100,'("    potential : ",F8.4)')wann_v(itype)
    write(100,*)
  endif
enddo
if (mpi_grid_root()) then
  if (wannier_lc) then
    write(100,*)
    write(100,'("number of linear combinations of WF : ",I4)')nwanlc
    do n=1,nwanlc
      write(100,'("  wf : ",I4)')n
      do i=1,wanlc_norb(n)
        write(100,'("    ",4I4)')(wanlc_iorb(j,i,n),j=1,4)
      enddo
    enddo  
  endif
  close(100)
endif
if (wannier_lc) nwantot=nwanlc

if (allocated(wanpos)) deallocate(wanpos)
allocate(wanpos(3,nwantot))
do j=1,nwantot
  if (wannier_lc) then
    write(*,'("TODO: positions for linear combinations of WFs")')
    call pstop
  else
    iatom=wan_info(1,j)
    wanpos(:,j)=atposc(:,ias2ia(iatom),ias2is(iatom))
  endif
enddo

if (allocated(wann_c)) deallocate(wann_c)
allocate(wann_c(nwantot,nstsv,nkptloc))
wann_c=zzero

! TODO: distribute over k-points
if (allocated(wann_h)) deallocate(wann_h)
allocate(wann_h(nwantot,nwantot,nkpt))
wann_h=zzero
if (allocated(wann_e)) deallocate(wann_e)
allocate(wann_e(nwantot,nkpt))
wann_e=0.d0
if (allocated(wann_p)) deallocate(wann_p)
allocate(wann_p(3,nwantot,nwantot,nkpt))
wann_p=zzero

if (allocated(wann_ene)) deallocate(wann_ene)
allocate(wann_ene(nwantot))
wann_ene=0.d0
if (allocated(wann_occ)) deallocate(wann_occ)
allocate(wann_occ(nwantot))
wann_occ=0.d0

if (allocated(wann_err_k)) deallocate(wann_err_k)
allocate(wann_err_k(nkptloc))
wann_err_k=0

return
end
