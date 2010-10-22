subroutine sic_init
use modmain
use mod_lf
implicit none
integer i,j,ias,n
integer i1,i2,i3,j1,j2,j3,is,ia,nt(3),ir0
real(8) d,vr0c,r0
logical l1
real(8) v1(3),v2(3),v3(3)
logical exist
logical, external :: vrinmt

call getnghbr(-0.1d0,wann_r_cutoff)
if (allocated(vtl)) deallocate(vtl)
allocate(vtl(3,maxvtl))
vtl=-1000000
ntr=0
tlim=0

do n=1,nwann
  ias=iwann(1,n)
  do j=1,nnghbr(ias)
    d=inghbr(2,j,ias)/1000000.d0
    if (d.le.wann_r_cutoff) then
      do i=1,3
        tlim(1,i)=min(tlim(1,i),inghbr(2+i,j,ias))
        tlim(2,i)=max(tlim(2,i),inghbr(2+i,j,ias))
      enddo
      l1=.true.
      do i=1,ntr
        if (all(vtl(:,i).eq.inghbr(3:5,j,ias))) l1=.false.
      enddo !i
      if (l1) then
        ntr=ntr+1
        if (ntr.gt.maxvtl) then
          write(*,'("Error(sic_init) : maxvtl is too small")')
          call pstop
        endif
        vtl(:,ntr)=inghbr(3:5,j,ias)
      endif
    endif
  enddo !j
enddo !n
do n=1,nwann
  ias=iwann(1,n)
  do i1=tlim(1,1)-1,tlim(2,1)+1
    do i2=tlim(1,2)-1,tlim(2,2)+1
      do i3=tlim(1,3)-1,tlim(2,3)+1
        v1(:)=i1*avec(:,1)+i2*avec(:,2)+i3*avec(:,3)-&
          atposc(:,ias2ia(ias),ias2is(ias))
        do j3=0,ngrid(3)-1
          v2(3)=dble(j3)/dble(ngrid(3))
          do j2=0,ngrid(2)-1
            v2(2)=dble(j2)/dble(ngrid(2))
            do j1=0,ngrid(1)-1
              v2(1)=dble(j1)/dble(ngrid(1))
              call r3mv(avec,v2,v3)
              v3(:)=v3(:)+v1(:)
              if (sqrt(sum(v3(:)**2)).le.wann_r_cutoff.and..not.&
                  vrinmt(v3,is,ia,nt,vr0c,ir0,r0)) then
                l1=.true.
                do i=1,ntr
                  if (all(vtl(:,i).eq.(/i1,i2,i3/))) l1=.false.
                enddo !i
                if (l1) then
                  ntr=ntr+1
                  if (ntr.gt.maxvtl) then
                    write(*,'("Error(sic_init) : maxvtl is too small")')
                    call pstop
                  endif
                  vtl(:,ntr)=(/i1,i2,i3/)
                endif
              endif
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
enddo
do i=1,3
  tlim(1,i)=minval(vtl(i,1:ntr))
  tlim(2,i)=maxval(vtl(i,1:ntr))
enddo
if (allocated(ivtit)) deallocate(ivtit)
allocate(ivtit(tlim(1,1):tlim(2,1),tlim(1,2):tlim(2,2),tlim(1,3):tlim(2,3)))
ivtit=-1
do i=1,ntr
  ivtit(vtl(1,i),vtl(2,i),vtl(3,i))=i
enddo
dim_t=dim2
ntrloc=mpi_grid_map(ntr,dim_t)
nwannloc=mpi_grid_map(nwann,dim_k)
if (mpi_grid_root()) then
  write(*,*)
  write(*,'("[sic_init] total number of translations : ",I3)')ntr
  write(*,'("[sic_init] local number of translations : ",I3)')ntrloc
  write(*,'("[sic_init] size of Wannier function arrays : ",I6," Mb")') &
    int(2*16.d0*(lmmaxvr*nrmtmax*natmtot+ngrtot)*ntrloc*nspinor*nwannloc/1048576.d0)
endif
call mpi_grid_barrier()
if (allocated(wvmt)) deallocate(wvmt)
allocate(wvmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor,nwannloc))
wvmt=zzero
if (allocated(wvir)) deallocate(wvir)  
allocate(wvir(ngrtot,ntrloc,nspinor,nwannloc))
wvir=zzero
if (allocated(wanmt)) deallocate(wanmt)
allocate(wanmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor,nwannloc))
wanmt=zzero
if (allocated(wanir)) deallocate(wanir)
allocate(wanir(ngrtot,ntrloc,nspinor,nwannloc))
wanir=zzero
if (allocated(sic_wann_e0)) deallocate(sic_wann_e0)
allocate(sic_wann_e0(nwann))
sic_wann_e0=0.d0
inquire(file="SIC_WANN_E0.OUT",exist=exist)
if (exist) then
  open(170,file="SIC_WANN_E0.OUT",form="FORMATTED",status="OLD")
  do n=1,nwann
    read(170,*)sic_wann_e0(n)
  enddo
  close(170)
endif
if (allocated(sic_wb)) deallocate(sic_wb)
allocate(sic_wb(nwann,nstfv,nspinor,nkptloc))
if (allocated(sic_wvb)) deallocate(sic_wvb)
allocate(sic_wvb(nwann,nstfv,nspinor,nkptloc))
if (allocated(sic_wann_h0k)) deallocate(sic_wann_h0k)
allocate(sic_wann_h0k(nwann,nwann,nkptloc))
sic_wann_h0k=zzero
sic_etot_correction=0.d0
tevecsv=.true.
return
end