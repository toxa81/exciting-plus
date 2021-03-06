subroutine sic_genpot(n,wanlm,wvlm,wanprop)
use modmain
use modxcifc
use mod_sic
use mod_hdf5
use mod_util
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: wanlm(lmmaxwan,s_nr,nspinor)
complex(8), intent(out) :: wvlm(lmmaxwan,s_nr,nspinor)
real(8), intent(out) :: wanprop(nwanprop)
! local variables
integer jr,ir,l,lm,ispn,lm1,lm2,lm3,lmmaxwanloc,lmloc
integer irloc,nrloc
real(8) t1
complex(8) zt1
real(8), allocatable :: rhotp(:,:,:)
real(8), allocatable :: rholm(:,:,:)
real(8), allocatable :: totrholm(:,:)
real(8), allocatable :: vhalm(:,:)
real(8), allocatable :: extp(:,:)
real(8), allocatable :: ectp(:,:)
real(8), allocatable :: vxtp(:,:,:)
real(8), allocatable :: vctp(:,:,:)
real(8), allocatable :: exclm(:,:)
real(8), allocatable :: vxclm(:,:,:)
real(8), allocatable :: rholm_tmp(:,:,:)
real(8), allocatable :: fx(:,:)
real(8), allocatable :: f1(:)
complex(8), allocatable :: wanlm_tmp(:,:,:)
real(8), allocatable :: vlm_tmp(:,:,:)
complex(8), allocatable :: wvlm_tmp(:,:,:)
complex(8), allocatable :: wantp(:,:)
real(8), external :: ddot
complex(8), external :: gauntyry
real(8), allocatable :: flm(:,:)
complex(8), allocatable :: zflm(:,:)
!
! TODO: generalize for non-collinear case; vxc will become 2x2 matrix
!

lmmaxwanloc=mpi_grid_map(lmmaxwan,dim_k)
allocate(f1(s_nr_min))

!
! 1) generate Rlm expansion of density
!
! wf expansion : w(r)=\sum_{L} w_{L}(r) Y_{L}(t,p)
! charge density : rho(r) = \sum_{L2} rho_{L2}(r) R_{L2}(t,p) = 
!   = \sum_{L1,L3} w_{L1}^{*}(r) Y_{L1}^{*}(t,p) * w_{L3}(r) Y_{L3} (t,p)
! rho_{L2}(r) = \sum_{L1,L3} w_{L1}^{*}(r) w_{L3}(r)   <Y_{L1} |R_{L2}| Y_{L3}>
allocate(wanlm_tmp(s_nr_min,lmmaxwan,nspinor))
allocate(rholm_tmp(s_nr_min,lmmaxwan,nspinor))
! rearrange wanlm in memory
do ispn=1,nspinor
  do lm=1,lmmaxwan
    wanlm_tmp(:,lm,ispn)=wanlm(lm,1:s_nr_min,ispn)
  enddo
enddo
rholm_tmp=0.d0
do lmloc=1,lmmaxwanloc
  lm2=mpi_grid_map(lmmaxwan,dim_k,loc=lmloc)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(flm,lm3,zt1,ispn)
  allocate(flm(s_nr_min,nspinor))
  flm(:,:)=0.d0
!$OMP DO
  do lm1=1,lmmaxwan
    do lm3=1,lmmaxwan
      zt1=gauntyry(lm2l(lm1),lm2l(lm2),lm2l(lm3),&
                  &lm2m(lm1),lm2m(lm2),lm2m(lm3))
      if (abs(zt1).gt.1d-12) then
        do ispn=1,nspinor
          flm(:,ispn)=flm(:,ispn)+dreal(dconjg(wanlm_tmp(:,lm1,ispn))*wanlm_tmp(:,lm3,ispn)*zt1)
        enddo !ispn
      endif
    enddo
  enddo
!$OMP END DO
!$OMP CRITICAL
  rholm_tmp(:,lm2,:)=rholm_tmp(:,lm2,:) + flm(:,:)
!$OMP END CRITICAL
  deallocate(flm)
!$OMP END PARALLEL
enddo
call mpi_grid_reduce(rholm_tmp(1,1,1),s_nr_min*lmmaxwan*nspinor,dims=(/dim_k/),all=.true.)
!deallocate(wanlm_tmp)
! total charge density of WF 
allocate(totrholm(s_nr_min,lmmaxwan))
totrholm=0.d0
do ispn=1,nspinor
  totrholm(:,:)=totrholm(:,:)+rholm_tmp(:,:,ispn)
enddo
! norm of rho
t1=rintegrate(s_nr_min,s_r,totrholm)
wanprop(wp_normrho)=t1*fourpi*y00
!
! 2) estimate the quadratic spread <r^2>-<r>^2
!
! Ry=-\frac{1}{2} \sqrt{\frac{3}{\pi }} \sin (t) \sin (p)  => y = -Ry * 2 * sqrt(pi/3)
! Rz= \frac{1}{2} \sqrt{\frac{3}{\pi }} \cos (t)           => z =  Rz * 2 * sqrt(pi/3)
! Rx=-\frac{1}{2} \sqrt{\frac{3}{\pi }} \sin (t) \cos (p)  => x = -Rx * 2 * sqrt(pi/3)
allocate(fx(s_nr_min,4))
do ir=1,s_nr_min
  fx(ir,1)=-2.d0*s_r(ir)*sqrt(pi/3.d0)*totrholm(ir,4)
  fx(ir,2)=-2.d0*s_r(ir)*sqrt(pi/3.d0)*totrholm(ir,2)
  fx(ir,3)=2.d0*s_r(ir)*sqrt(pi/3.d0)*totrholm(ir,3)
  fx(ir,4)=(s_r(ir)**2)*totrholm(ir,1)/y00
enddo
wanprop(wp_spread_x)=rintegrate(s_nr_min,s_r,fx(1,1))
wanprop(wp_spread_y)=rintegrate(s_nr_min,s_r,fx(1,2))
wanprop(wp_spread_z)=rintegrate(s_nr_min,s_r,fx(1,3))
wanprop(wp_spread)=rintegrate(s_nr_min,s_r,fx(1,4))-&
  &(wanprop(wp_spread_x)**2+wanprop(wp_spread_y)**2+wanprop(wp_spread_z)**2)
deallocate(fx)
!
! 3) compute Hartree potential
!
allocate(vhalm(lmmaxwan,s_nr_min))
vhalm=0.d0
do lmloc=1,lmmaxwanloc
  lm=mpi_grid_map(lmmaxwan,dim_k,loc=lmloc)
  l=lm2l(lm)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jr,t1)
  do ir=1,s_nr_min
    t1=0.d0
    do jr=1,ir
      t1=t1+((s_r(jr)/s_r(ir))**l)*totrholm(jr,lm)*s_rw(jr)/s_r(ir)
    enddo
    do jr=ir+1,s_nr_min
      t1=t1+((s_r(ir)/s_r(jr))**l)*totrholm(jr,lm)*s_rw(jr)/s_r(jr)
    enddo
    vhalm(lm,ir)=t1*fourpi/(2*l+1)
  enddo
!$OMP END PARALLEL DO
enddo
call mpi_grid_reduce(vhalm(1,1),lmmaxwan*s_nr_min,dims=(/dim_k/),all=.true.)
deallocate(totrholm)
allocate(rholm(lmmaxwan,s_nr_min,nspinor))
do ispn=1,nspinor
  do lm=1,lmmaxwan
    rholm(lm,:,ispn)=rholm_tmp(:,lm,ispn)
  enddo
enddo
deallocate(rholm_tmp)
! compute vha=<V_h|rho>
f1=0.d0
do ir=1,s_nr_min
  do ispn=1,nspinor
    f1(ir)=f1(ir)+ddot(lmmaxwan,rholm(1,ir,ispn),1,vhalm(1,ir),1)
  enddo
enddo
wanprop(wp_vha)=rintegrate(s_nr_min,s_r,f1)
!
! 4) compute XC potential
!
! compute density 
allocate(rhotp(s_ntp,s_nr_min,nspinor))
allocate(wantp(s_ntp,s_nr_min))
do ispn=1,nspinor
  call zgemm('T','N',s_ntp,s_nr_min,lmmaxwan,zone,s_ylmf,lmmaxwan,&
    &wanlm(1,1,ispn),lmmaxwan,zzero,wantp,s_ntp)
  rhotp(:,1:s_nr_min,ispn)=abs(wantp(:,1:s_nr_min))**2
enddo
deallocate(wantp)

nrloc=mpi_grid_map(s_nr_min,dim_k)

allocate(extp(s_ntp,s_nr_min))
extp=0.d0
allocate(ectp(s_ntp,s_nr_min))
ectp=0.d0
allocate(vxtp(s_ntp,s_nr_min,nspinor))
vxtp=0.d0
allocate(vctp(s_ntp,s_nr_min,nspinor))
vctp=0.d0

do irloc=1,nrloc
  ir=mpi_grid_map(s_nr_min,dim_k,loc=irloc)
  if (spinpol) then
    call xcifc(xctype,n=s_ntp,rhoup=rhotp(1,ir,1),rhodn=rhotp(1,ir,2),&
      &ex=extp(1,ir),ec=ectp(1,ir),vxup=vxtp(1,ir,1),vxdn=vxtp(1,ir,2),&
      &vcup=vctp(1,ir,1),vcdn=vctp(1,ir,2))
  else
    call xcifc(xctype,n=s_ntp,rho=rhotp(1,ir,1),ex=extp(1,ir),ec=ectp(1,ir),&
      &vx=vxtp(1,ir,1),vc=vctp(1,ir,1))
  endif
enddo
call mpi_grid_reduce(extp(1,1),s_ntp*s_nr_min,dims=(/dim_k/),all=.true.)
call mpi_grid_reduce(ectp(1,1),s_ntp*s_nr_min,dims=(/dim_k/),all=.true.)
call mpi_grid_reduce(vxtp(1,1,1),s_ntp*s_nr_min*nspinor,dims=(/dim_k/),all=.true.)
call mpi_grid_reduce(vctp(1,1,1),s_ntp*s_nr_min*nspinor,dims=(/dim_k/),all=.true.)

deallocate(rhotp)

allocate(exclm(lmmaxwan,s_nr_min))
! compute exchange and correlation contributions to the XC energy
call sic_rbsht(s_nr_min,extp,exclm)
f1=0.d0
do ir=1,s_nr_min
  do ispn=1,nspinor
    f1(ir)=f1(ir)+ddot(lmmaxwan,exclm(1,ir),1,rholm(1,ir,ispn),1)
  enddo
enddo
wanprop(wp_ex)=rintegrate(s_nr_min,s_r,f1)
call sic_rbsht(s_nr_min,ectp,exclm)
f1=0.d0
do ir=1,s_nr_min
  do ispn=1,nspinor
    f1(ir)=f1(ir)+ddot(lmmaxwan,exclm(1,ir),1,rholm(1,ir,ispn),1)
  enddo
enddo
wanprop(wp_ec)=rintegrate(s_nr_min,s_r,f1)
if (sicec) then
  wanprop(wp_exc)=wanprop(wp_ex)+wanprop(wp_ec)
else
  wanprop(wp_exc)=wanprop(wp_ex)
endif
deallocate(extp,ectp,exclm)
allocate(vxclm(lmmaxwan,s_nr_min,nspinor))
! save XC potential in vxtp and expand in real spherical harmonics   
do ispn=1,nspinor
  if (sicvc) vxtp(:,:,ispn)=vxtp(:,:,ispn)+vctp(:,:,ispn)
  call sic_rbsht(s_nr_min,vxtp(1,1,ispn),vxclm(1,1,ispn))
enddo
!! write XC potential
!if (sic_write_vxclm.and.mpi_grid_root()) then
!  write(path,'("/wannier_functions/",I4.4)')n
!  call hdf5_write("sic.hdf5",path,"vxclm",vxclm(1,1,1),&
!    &(/lmmaxwan,s_nr_min,nspinor/))
!endif
! compute vxc=<W_n|V_xc|W_n>; in the collinear case this is 
!  \sum_{\sigma} <V_xc^{\sigma}|rho^{\sigma}>
f1=0.d0
do ir=1,s_nr_min
  do ispn=1,nspinor
    f1(ir)=f1(ir)+ddot(lmmaxwan,rholm(1,ir,ispn),1,vxclm(1,ir,ispn),1)
  enddo
enddo
wanprop(wp_vxc)=rintegrate(s_nr_min,s_r,f1)
deallocate(rholm,vxtp,vctp)

! compute <V_n|rho>
wanprop(wp_vsic)=wanprop(wp_vha)+wanprop(wp_vxc)
! add Hartree potential to XC
do ispn=1,nspinor
  vxclm(:,:,ispn)=vxclm(:,:,ispn)+vhalm(:,:)
enddo
deallocate(vhalm)
!! write total SIC potential
!if (sic_write_vlm.and.mpi_grid_root()) then
!  write(path,'("/wannier_functions/",I4.4)')n
!  call hdf5_write("sic.hdf5",path,"vlm",vxclm(1,1,1),&
!    &(/lmmaxwan,s_nr_min,nspinor/))
!endif

!
! 5) multiply Wannier function with potential and change sign
!
call timer_start(t_sic_wvprod)
allocate(vlm_tmp(s_nr_min,lmmaxwan,nspinor))
allocate(wvlm_tmp(s_nr_min,lmmaxwan,nspinor))
! rearrange array in memory
do ispn=1,nspinor
  do lm=1,lmmaxwan
    vlm_tmp(:,lm,ispn)=vxclm(lm,1:s_nr_min,ispn)
  enddo
enddo
wvlm_tmp=zzero
do lmloc=1,lmmaxwanloc
  lm3=mpi_grid_map(lmmaxwan,dim_k,loc=lmloc)
!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(zflm,lm2,zt1,ispn)
  allocate(zflm(s_nr_min,nspinor))
  zflm=zzero
!$OMP DO  
  do lm1=1,lmmaxwan
    do lm2=1,lmmaxwan
      zt1=gauntyry(lm2l(lm3),lm2l(lm2),lm2l(lm1),&
                  &lm2m(lm3),lm2m(lm2),lm2m(lm1))
      if (abs(zt1).gt.1d-12) then
        do ispn=1,nspinor
          zflm(:,ispn)=zflm(:,ispn)-wanlm_tmp(:,lm1,ispn)*vlm_tmp(:,lm2,ispn)*zt1
        enddo !ispn
      endif
    enddo
  enddo
!$OMP END DO
!$OMP CRITICAL
  wvlm_tmp(:,lm3,:)=wvlm_tmp(:,lm3,:) + zflm(:,:)
!$OMP END CRITICAL
  deallocate(zflm)
!$OMP END PARALLEL
enddo
call mpi_grid_reduce(wvlm_tmp(1,1,1),s_nr_min*lmmaxwan*nspinor,dims=(/dim_k/),all=.true.) 
wvlm=zzero
do ispn=1,nspinor
  do lm=1,lmmaxwan
    wvlm(lm,1:s_nr_min,ispn)=wvlm_tmp(1:s_nr_min,lm,ispn)
  enddo
enddo
deallocate(vlm_tmp,wvlm_tmp,wanlm_tmp,vxclm,f1)
call timer_stop(t_sic_wvprod)
end subroutine   

