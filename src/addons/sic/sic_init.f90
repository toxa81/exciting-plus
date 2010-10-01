subroutine sic_init
use modmain
use mod_lf
implicit none
integer i,j,ias,n
integer i1,i2,i3,j1,j2,j3
real(8) d
logical l1
real(8) v1(3),v2(3),v3(3)

call getnghbr(-0.1d0,wann_r_cutoff)
if (allocated(wann_r_nnghbr)) deallocate(wann_r_nnghbr)
allocate(wann_r_nnghbr(natmtot))
wann_r_nnghbr=nnghbr
n=maxval(nnghbr)
if (allocated(wann_r_inghbr)) deallocate(wann_r_inghbr)
allocate(wann_r_inghbr(5,n,natmtot))
do i=1,natmtot
  wann_r_inghbr(1:5,1:nnghbr(i),i)=inghbr(1:5,1:nnghbr(i),i)
enddo

tlim=0
do n=1,nwann
  ias=iwann(1,n)
  do j=1,wann_r_nnghbr(ias)
    d=wann_r_inghbr(2,j,ias)/1000000.d0
    if (d.le.wann_r_cutoff) then
      do i=1,3
        tlim(1,i)=min(tlim(1,i),wann_r_inghbr(2+i,j,ias))
        tlim(2,i)=max(tlim(2,i),wann_r_inghbr(2+i,j,ias))
      enddo
    endif
  enddo !j
enddo !n
if (allocated(vtl)) deallocate(vtl)
allocate(vtl(3,maxvtl))
vtl=-10000
ntr=0
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
              if (sqrt(sum(v3(:)**2)).le.wann_r_cutoff) then
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
sic_etot_correction=0.d0
if (allocated(wvmt)) deallocate(wvmt)
allocate(wvmt(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor,nwann))
if (allocated(wvir)) deallocate(wvir)  
allocate(wvir(ngrtot,ntrloc,nspinor,nwann))
if (allocated(hmltsv)) deallocate(hmltsv)
allocate(hmltsv(nstsv,nstsv,nkptloc))
return
end