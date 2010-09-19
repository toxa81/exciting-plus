subroutine genmegq(iq,tout,tg0q)
use modmain
use mod_nrkp
use mod_addons_q
implicit none
integer, intent(in) :: iq
logical, intent(in) :: tout
logical, intent(in) :: tg0q

! allocatable arrays
integer, allocatable :: igkignr_jk(:)
complex(8), allocatable :: wfsvmt_jk(:,:,:,:,:)
complex(8), allocatable :: wfsvit_jk(:,:,:)
integer ngknr_jk
integer i,ikstep,sz,ig
integer nkstep,ik,ist1,ist2,ikloc
real(8) t1,t2,t3,t4,t5,dn1
integer lmaxexp,lmmaxexp
integer np
character*100 :: qnm,qdir,fout

! comment:
! the subroutine computes <psi_{n,k}|e^{-i(G+q)x}|psi_{n',k+q}>  and
!  <W_n|e^{-i(G+q)x}|W_{n'T}> 

! maximum l for exponent expansion
lmaxexp=lmaxvr
lmmaxexp=(lmaxexp+1)**2
call getqdir(iq,vqm(:,iq),qdir)
call getqname(vqm(:,iq),qnm)
qnm=trim(qdir)//"/"//trim(qnm)
wproc=.false.
if (mpi_grid_root((/dim_k,dim_b/)).and.tout) then
  wproc=.true.
  fout=trim(qnm)//"_ME.OUT"
  open(150,file=trim(fout),form='formatted',status='replace')
endif

if (wproc) then
  write(150,*)
  write(150,'("Calculation of matrix elements:")')
  write(150,'("  <n,k|e^{-i(G+q)x}|n'',k+q>")')
endif
call timer_start(1,reset=.true.)
! initialize G+q vector arays
call init_gq(iq,lmaxexp,lmmaxexp,tg0q)
! initialize k+q array
call init_kq(iq)
! initialize interband transitions
call init_band_trans
! initialize Gaunt-like coefficients 
call init_gntuju(iq,lmaxexp)
call timer_stop(1)
if (wproc) then
  write(150,*)
  write(150,'("maximum |G+q| [1/a.u.]                        : ",G18.10)')gqmax  
  write(150,'("number of G-vectors                           : ",I4)')ngvecme   
  write(150,*)
  write(150,'("q-vector (lat.coord.)                         : ",&
    & 3G18.10)')vqlnr(:,iq)
  write(150,'("q-vector (Cart.coord.) [1/a.u.]               : ",&
    & 3G18.10)')vqcnr(:,iq)
  t1=sqrt(vqcnr(1,iq)**2+vqcnr(2,iq)**2+vqcnr(3,iq)**2)
  write(150,'("q-vector length [1/a.u.]                      : ",&
    & G18.10)')t1
  write(150,'("q-vector length [1/A]                         : ",&
    & G18.10)')t1/au2ang
  write(150,'("G-vector to reduce q to first BZ (lat.coord.) : ",&
    & 3I4)')ivg(:,ig0q(iq))
  write(150,'("index of G-vector                             : ",&
    & I4)')ig0q(iq)
  write(150,'("reduced q-vector (lat.coord.)                 : ",&
    & 3G18.10)')vql(:,iq)
  write(150,'("reduced q-vector (Cart.coord.) [1/a.u.]       : ",&
    & 3G18.10)')vqc(:,iq)
  write(150,*)
  write(150,'("Bloch functions band interval (N1,N2 or E1,E2) : ",2F8.3)')&
    chi0_include_bands(1),chi0_include_bands(2)
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
allocate(wfsvmt_jk(lmmaxvr,nufrmax,natmtot,nspinor,nstsv))
allocate(wfsvit_jk(ngkmax,nspinor,nstsv))
allocate(igkignr_jk(ngkmax))
if (wannier_megq) then
  if (allocated(megqwan)) deallocate(megqwan)
  allocate(megqwan(nmegqwan,ngvecme))
  megqwan(:,:)=zzero
  if (allocated(wann_c_jk)) deallocate(wann_c_jk)
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
  call getwfkq(ikstep,ngknr_jk,igkignr_jk,wfsvmt_jk,wfsvit_jk)
  call timer_stop(1)
! compute matrix elements  
  call timer_start(2)
  if (ikstep.le.nkptnrloc) then
    call genmegqblh(iq,ikstep,ngknr(ikstep),ngknr_jk,igkignr(1,ikstep),&
      igkignr_jk,wfsvmtnrloc(1,1,1,1,1,ikstep),wfsvmt_jk,&
      wfsvitnrloc(1,1,1,ikstep),wfsvit_jk)
  endif !ikstep.le.nkptnrloc
  call timer_stop(2)
enddo !ikstep
!call printmegqblh(iq)
! for G=q=0: e^{iqx}=1+iqx
! from "Formalism of Bnad Theory" by E.I. Blount:
!   v=p/m
!   v=-i/\hbar [x,H]
! -i(xH-Hx)=p 
! xH-Hx=ip
! <nk|xH|n'k>-<nk|Hx|n'k>=E_{n'k}<nk|x|n'k>-E_{nk}<nk|x|n'k>
! <nk|x|n'k>*(E_{n'k}-E_{nk})=i<nk|p|n'k>
! <nk|x|n'k>=i<nk|p|n'k>/(E_{n'k}-E_{nk})
! <nk|e^{iqx}|n'k>=<nk|1+iqx|n'k>=\delta_{nn'}+iq<nk|x|n'k>
! <nk|e^{iqx}|n'k>=\delta_{nn'}-q*<nk|p|n'k>/(E_{n'k}-E_{nk})
if (all(vqm(:,iq).eq.0).and.allocated(pmatnrloc)) then
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do ig=1,ngvecme
      if (igqig(ig,iq).eq.1) then
        megqblh(:,ig,ikloc)=zzero
        do i=1,nmegqblhloc(1,ikloc)
          ist1=bmegqblh(1,i+nmegqblhloc(2,ikloc),ikloc)
          ist2=bmegqblh(2,i+nmegqblhloc(2,ikloc),ikloc)
          t1=evalsvnr(ist2,ik)-evalsvnr(ist1,ik)
          if (ist1.eq.ist2) megqblh(i,ig,ikloc)=zone
          if (abs(t1).gt.1d-8) then
            megqblh(i,ig,ikloc)=megqblh(i,ig,ikloc)-&
              dot_product(vq0c(:,iq),pmatnrloc(:,ist1,ist2,ikloc))/t1          
          endif
        enddo
      endif
    enddo
  enddo
endif
!call printmegqblh(iq+100)
if (wannier_megq) then
  call timer_start(6,reset=.true.)
! compute matrix elements of e^{-i(G+q)x} in the basis of Wannier functions
  call genmegqwan(iq)
! sum over all k-points and interband transitions to get <n,T=0|e^{-i(G+q)x}|n',T'>
  call mpi_grid_reduce(megqwan(1,1),nmegqwan*ngvecme,dims=(/dim_k,dim_b/),&
    all=.true.)
  megqwan=megqwan/nkptnr
  call timer_stop(6)
  if (wproc) then
    write(150,*)
    write(150,'("Time for megqwan : ",F8.2)')timer_get_value(6)
  endif
  call printmegqwan(iq)
endif
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
deallocate(wfsvmt_jk)
deallocate(wfsvit_jk)
deallocate(igkignr_jk)
deallocate(ngntuju)
deallocate(igntuju)
deallocate(gntuju)
call mpi_grid_barrier((/dim_k,dim_b/))
if (wproc) then
  write(150,*)
  write(150,'("Done.")')
  call flushifc(150)
endif
return
end
