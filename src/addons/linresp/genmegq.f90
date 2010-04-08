#ifdef _HDF5_
subroutine genmegq(ivq0m)
use modmain
use mod_nrkp
implicit none
! q-vector in k-mesh coordinates
integer, intent(in) :: ivq0m(3)

! allocatable arrays
integer, allocatable :: igkignr_jk(:)
complex(8), allocatable :: wfsvmt_jk(:,:,:,:,:)
complex(8), allocatable :: wfsvit_jk(:,:,:)
complex(8), allocatable :: wann_c_jk(:,:,:)
integer ngknr_jk
integer i,ikstep,sz,complete
integer nkstep
real(8) t1,t2,t3,t4,t5,dn1
integer lmaxexp,lmmaxexp
integer np
character*100 :: qnm,fout,fme,fu
logical exist

! comment:
! the subroutine computes <psi_{n,k}|e^{-i(G+q)x}|psi_{n',k+q}> 
! 
! switch write_megq_file controls the reading and writing of ME file
! when we write ME we have two choices: write to single file or write
!  to multiple files

! maximum l for exponent expansion
lmaxexp=lmaxvr+2
lmmaxexp=(lmaxexp+1)**2

call qname(ivq0m,qnm)
qnm="./qv/"//trim(qnm)//"/"//trim(qnm)
wproc=.false.
if (mpi_grid_root((/dim_k,dim_b/))) then
  wproc=.true.
  fout=trim(qnm)//"_ME.OUT"
  open(150,file=trim(fout),form='formatted',status='replace')
endif

complete=0
fme=trim(qnm)//"_me.hdf5"
if (mpi_grid_root((/dim_k,dim_b/))) then
  inquire(file=trim(fme),exist=exist)
  if (exist) then
    call read_integer(complete,1,trim(fme),'/parameters','complete')
  endif
endif
call mpi_grid_bcast(complete,dims=(/dim_k,dim_b/))
if (complete.eq.1) goto 30

if (crpa) then
  if (mpi_grid_root((/dim_k,dim2/))) then
    fu=trim(qnm)//"_U"
    inquire(file=trim(fu),exist=exist)
  endif
  call mpi_grid_bcast(exist,dims=(/dim_k,dim2/))
  if (exist) goto 30
endif

if (wproc) then
  write(150,*)
  write(150,'("Calculation of matrix elements:")')
  write(150,'("  <n,k|e^{-i(G+q)x}|n'',k+q>")')
endif

call timer_start(1,reset=.true.)
! initialize G, q and G+q vectors
call init_g_q_gq(ivq0m,lmaxexp,lmmaxexp)
! initialize k+q array
call init_kq
! initialize interband transitions
call init_band_trans
! initialize Gaunt-like coefficients 
call init_gntuju(lmaxexp)
call timer_stop(1)

if (wproc) then
  write(150,*)
  write(150,'("G-shell limits      : ",2I4)')gshme1,gshme2
  write(150,'("G-vector limits     : ",2I4)')gvecme1,gvecme2
  write(150,'("number of G-vectors : ",I4)')ngvecme   
  write(150,*)
  write(150,'("q-vector (lat.coord.)                        : ",&
    & 3G18.10)')vq0l
  write(150,'("q-vector (Cart.coord.) [a.u.]                : ",&
    & 3G18.10)')vq0c
  write(150,'("q-vector length [a.u.]                       : ",&
    & G18.10)')sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)
  write(150,'("q-vector length [1/A]                        : ",&
    & G18.10)')sqrt(vq0c(1)**2+vq0c(2)**2+vq0c(3)**2)/au2ang
  write(150,'("G-vector to reduce q to first BZ (lat.coord.): ",&
    & 3I4)')ivg(:,lr_igq0)
  write(150,'("index of G-vector                            : ",&
    & I4)')lr_igq0
  write(150,'("reduced q-vector (lat.coord.)                : ",&
    & 3G18.10)')vq0rl
  write(150,'("reduced q-vector (Cart.coord.) [a.u.]        : ",&
    & 3G18.10)')vq0rc
  write(150,*)
  write(150,'("Bloch functions band interval (N1,N2 or E1,E2) : ",2F8.3)')&
    lr_e1,lr_e2
  write(150,*)
  write(150,'("Minimal energy transition (eV) : ",F12.6)')lr_min_e12*ha2ev    
  write(150,*)
  write(150,'("Maximum number of interband transitions : ",I5)')nmegqblhmax
  if (wannier_megq) then
    write(150,*)
    write(150,'("Maximum number of interband transitions for megqwan : ",I5)')nmegqblhwanmax
  endif
  sz=int(16.d0*ngvecme*nmegqblhlocmax*nkptnrloc/1048576.d0)
  write(150,*)
  write(150,'("Array size of matrix elements in Bloch basis (MB) : ",I6)')sz
  if (wannier_megq) then
    sz=int(16.d0*nmegqwan*ngvecme/1048576.d0)
    write(150,*)
    write(150,'("Array size of matrix elements in Wannier basis (MB) : ",I6)')sz
  endif   
  sz=int(24.d0*ngntujumax*natmcls*ngvecme/1048576.d0)
  write(150,*)
  write(150,'("Maximum number of Gaunt-like coefficients : ",I8)')ngntujumax
  write(150,'("Array size of Gaunt-like coefficients (MB) : ",I6)')sz
  write(150,*)
  write(150,'("Init done in ",F8.2," seconds")')timer_get_value(1)
  call flushifc(150)
endif

if (allocated(megqblh)) deallocate(megqblh)
allocate(megqblh(nmegqblhlocmax,ngvecme,nkptnrloc))
megqblh(:,:,:)=zzero
if (wannier_megq) then
  if (allocated(megqwan)) deallocate(megqwan)
  allocate(megqwan(nmegqwan,ngvecme))
  megqwan(:,:)=zzero
endif

if (write_megq_file) call write_me_header(qnm)

allocate(wfsvmt_jk(lmmaxvr,nrfmax,natmtot,nspinor,nstsv))
allocate(wfsvit_jk(ngkmax,nspinor,nstsv))
allocate(igkignr_jk(ngkmax))
if (wannier_megq) then
  allocate(wann_c_jk(nwann,nstsv,nkptnrloc))
endif

i=0
nkstep=mpi_grid_map(nkptnr,dim_k,x=i)
call timer_reset(1)
call timer_reset(2)
call timer_reset(3)
call timer_reset(4)
call timer_reset(5)
do ikstep=1,nkstep
! transmit wave-functions
  call timer_start(1)
  call getwfkq(ikstep,ngknr_jk,igkignr_jk,wfsvmt_jk,wfsvit_jk,wann_c_jk)
  call timer_stop(1)
! compute matrix elements  
  call timer_start(2)
  if (ikstep.le.nkptnrloc) then
    call genmegqblh(ikstep,ngknr(ikstep),ngknr_jk,igkignr(1,ikstep),&
      igkignr_jk,wfsvmtloc(1,1,1,1,1,ikstep),wfsvmt_jk,&
      wfsvitloc(1,1,1,ikstep),wfsvit_jk)
! add contribution from k-point to the matrix elements of e^{-i(G+q)x} in 
!  the basis of Wannier functions
    if (wannier_megq) then
      call genmegqwan(ikstep,wann_c_jk)
    endif !wannier
  endif !ikstep.le.nkptnrloc
  call timer_stop(2)
enddo !ikstep
! time for wave-functions send/recieve
t1=timer_get_value(1)
call mpi_grid_reduce(t1,dims=(/dim_k,dim_b/))
! total time for matrix elements calculation
t2=timer_get_value(2)
call mpi_grid_reduce(t2,dims=(/dim_k,dim_b/))
! time to precompute MT
t3=timer_get_value(3)
call mpi_grid_reduce(t3,dims=(/dim_k,dim_b/))
! time to precompute IT
t4=timer_get_value(4)
call mpi_grid_reduce(t4,dims=(/dim_k,dim_b/))
! time to compute ME
t5=timer_get_value(5)
call mpi_grid_reduce(t5,dims=(/dim_k,dim_b/))
! approximate number of matrix elements
dn1=1.d0*nmegqblhmax*ngvecme*nkptnr
if (wannier_megq) dn1=dn1+1.d0*nmegqwan*ngvecme
np=mpi_grid_size(dim_k)*mpi_grid_size(dim_b)
if (wproc) then
  write(150,*)
  write(150,'("Average time (seconds/proc)")')
  write(150,'("  send and receive wave-functions  : ",F8.2)')t1/np
  write(150,'("  compute matrix elements          : ",F8.2)')t2/np
  write(150,'("    precompute muffin-tin part     : ",F8.2)')t3/np
  write(150,'("    precompute interstitial part   : ",F8.2)')t4/np
  write(150,'("    multiply wave-functions        : ",F8.2)')t5/np
  write(150,'("Speed (me/sec/proc)                : ",F10.2)')dn1/t2
  call flushifc(150)
endif

if (wannier_megq) then
! sum over all k-points and interband transitions to get <n,T=0|e^{-i(G+q)x|n',T'>
  call mpi_grid_reduce(megqwan(1,1),nmegqwan*ngvecme,dims=(/dim_k,dim_b/),&
    all=.true.)
  megqwan=megqwan/nkptnr
endif

if (write_megq_file) then
  if (wproc) then
    write(150,*)
    write(150,'("Writing matrix elements")')
  endif
  call timer_start(3,reset=.true.)
  call write_me(qnm)
  call timer_stop(3)
  if (wproc) write(150,'(" Done in : ",F8.2)')timer_get_value(3)
endif  

deallocate(wfsvmt_jk)
deallocate(wfsvit_jk)
deallocate(igkignr_jk)
if (wannier_megq) then
  deallocate(wann_c_jk)
endif
deallocate(ngntuju)
deallocate(igntuju)
deallocate(gntuju)

if (mpi_grid_root((/dim_k,dim_b/)).and.write_megq_file) then
  complete=1
  call rewrite_integer(complete,1,trim(fme),'/parameters','complete')
endif

call mpi_grid_barrier((/dim_k,dim_b/))

30 continue
if (wproc) then
  write(150,*)
  write(150,'("Done.")')
  call flushifc(150)
endif

return
end
#endif
