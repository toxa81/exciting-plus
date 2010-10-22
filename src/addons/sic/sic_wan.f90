subroutine sic_wan(fout)
use modmain
use mod_lf
use mod_nrkp
implicit none
! arguments
integer, intent(in) :: fout
! local variables
complex(8), allocatable :: wanmt0(:,:,:,:,:)
complex(8), allocatable :: wanir0(:,:,:)
complex(8), allocatable :: ovlp(:)
integer itloc,it,n,ispn,vl(3),nloc,n1loc,h,h1,n1,i,j
real(8) t1,t2,t3
complex(8) z1

if (wproc) then
  write(fout,*)
  write(fout,'("sic_wan.f90")')
  write(fout,'("generate Wannier functions on a grid")')
  write(fout,'(80("-"))')
endif
call timer_reset(1)
call timer_reset(2)
allocate(wanmt0(lmmaxvr,nrmtmax,natmtot,nspinor,nwannloc))
allocate(wanir0(ngrtot,nspinor,nwannloc))
do itloc=1,ntrloc
  it=mpi_grid_map(ntr,dim_t,loc=itloc)
! generate Wannier functions on a mesh
  call sic_genwann(vtl(1,it),ngknr,igkignr,wanmt0,wanir0)
  do ispn=1,nspinor
    do n=1,nwannloc
      wanmt(:,:,:,itloc,ispn,n)=wanmt0(:,:,:,ispn,n)
      wanir(:,itloc,ispn,n)=wanir0(:,ispn,n)
    enddo !n
  enddo !ispn
enddo !it
deallocate(wanmt0,wanir0)
deallocate(wann_unkmt)
deallocate(wann_unkit)
if (wproc) then
  write(fout,*)
  write(fout,'("time for muffin-tin part (sec.)   : ",F8.3)')timer_get_value(1)
  write(fout,'("time for interstitial part (sec.) : ",F8.3)')timer_get_value(2)
  call flushifc(fout)
endif

allocate(wanmt0(lmmaxvr,nrmtmax,natmtot,ntrloc,nspinor))
allocate(wanir0(ngrtot,ntrloc,nspinor))
allocate(ovlp(nmegqwan))
ovlp=zzero
! compute overlap integrals
do i=1,nmegqwan
  n=imegqwan(1,i)
  nloc=mpi_grid_map(nwann,dim_k,x=h,glob=n)
  n1=imegqwan(2,i)
  n1loc=mpi_grid_map(nwann,dim_k,x=h1,glob=n1)
  vl(:)=imegqwan(3:5,i)
  if (mpi_grid_x(dim_k).eq.h1) then
    call mpi_grid_send(wanmt(1,1,1,1,1,n1loc),lmmaxvr*nrmtmax*natmtot*ntrloc*nspinor,&
      dims=(/dim_k/),dest=(/h/),tag=1)
    call mpi_grid_send(wanir(1,1,1,n1loc),ngrtot*ntrloc*nspinor,&
      dims=(/dim_k/),dest=(/h/),tag=2)
  endif
  if (mpi_grid_x(dim_k).eq.h) then
    call mpi_grid_recieve(wanmt0(1,1,1,1,1),lmmaxvr*nrmtmax*natmtot*ntrloc*nspinor,&
      dims=(/dim_k/),src=(/h1/),tag=1)
    call mpi_grid_recieve(wanir0(1,1,1),ngrtot*ntrloc*nspinor,&
      dims=(/dim_k/),src=(/h1/),tag=2)
  endif
  if (mpi_grid_x(dim_k).eq.h) then
    do ispn=1,nspinor
      ovlp(i)=ovlp(i)+lf_dot_lf(.true.,wanmt(1,1,1,1,ispn,nloc),&
        wanir(1,1,ispn,nloc),vl,wanmt0(1,1,1,1,ispn),wanir0(1,1,ispn))
    enddo
  endif
enddo
call mpi_grid_reduce(ovlp(1),nmegqwan,dims=(/dim_k/))
! check orthonormality
t1=0.d0
t2=0.d0
t3=0.d0
do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  vl(:)=imegqwan(3:5,i)
  j=idxmegqwan(n1,n,-vl(1),-vl(2),-vl(3))
  z1=ovlp(i)
  if (n.eq.n1.and.vl(1).eq.0.and.vl(2).eq.0.and.vl(3).eq.0) then
    z1=z1-zone
  endif
  t2=max(t2,abs(z1))
  t1=t1+abs(z1)
  t3=max(t3,abs(ovlp(i)-dconjg(ovlp(j))))
enddo
if (wproc) then
  write(fout,*)
  write(fout,'("Wannier overlap integrals (n n1  T  <w_n|w_{n1,T}>)")')
  do i=1,nmegqwan
    write(151,'(I4,4X,I4,4X,3I3,4X,2G18.10)')imegqwan(:,i),&
      dreal(ovlp(i)),dimag(ovlp(i))
  enddo
  write(fout,*)
  write(fout,'("Maximum deviation from norm                 : ",F12.6)')t2
!  write(fout,'("Average deviation from norm : ",F12.6)')t1/nmegqwan
  write(fout,'("Maximum of <w_n|w_{n1,T}> - <w_n1|w_{n,-T}> : ",G18.10)')t3
  call timestamp(fout,"done with Wannier functions")
  call flushifc(151)
endif
deallocate(wanmt0,wanir0,ovlp)
return
end