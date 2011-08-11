subroutine sic_genpot(n,wanlm,wvlm,wanprop)
use modmain
use modxcifc
use mod_sic
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: wanlm(lmmaxwan,s_nr,nspinor)
complex(8), intent(out) :: wvlm(lmmaxwan,s_nr,nspinor)
real(8), intent(out) :: wanprop(nwanprop)
! local variables
integer jr,ir,l,lm,ispn,lm1,lm2,lm3,lmmaxwanloc,lmloc
real(8) t1,x(3),x2
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
real(8), allocatable :: f2(:,:,:)
complex(8), allocatable :: f1(:,:,:)
complex(8), allocatable :: f3(:,:,:)
complex(8), allocatable :: wantp(:,:)
real(8), external :: ddot
complex(8), external :: gauntyry
!
! TODO: generalize for non-collinear case; vxc will become 2x2 matrix
!
allocate(rhotp(s_ntp,s_nr_min,nspinor))
allocate(rholm(lmmaxwan,s_nr_min,nspinor))
allocate(totrholm(lmmaxwan,s_nr_min))
allocate(vhalm(lmmaxwan,s_nr_min))
allocate(extp(s_ntp,s_nr_min))
allocate(ectp(s_ntp,s_nr_min))
allocate(exclm(lmmaxwan,s_nr_min))
allocate(vxtp(s_ntp,s_nr_min,nspinor))
allocate(vctp(s_ntp,s_nr_min,nspinor))
allocate(vxclm(lmmaxwan,s_nr_min,nspinor))

lmmaxwanloc=mpi_grid_map2(lmmaxwan,dims=(/dim_k,dim2/))

! TODO: memory optimizaton?
totrholm=0.d0
rholm=0.d0
allocate(wantp(s_ntp,s_nr_min))
do ispn=1,nspinor
  call zgemm('T','N',s_ntp,s_nr_min,lmmaxwan,zone,s_ylmf,lmmaxwan,&
    wanlm(1,1,ispn),lmmaxwan,zzero,wantp,s_ntp)
  rhotp(:,1:s_nr_min,ispn)=abs(wantp(:,1:s_nr_min))**2
enddo
deallocate(wantp)
! charge density
! w(r)=\sum_{L} w_{L}(r) Y_{L}(t,p)
! rho(r) = \sum_{L2} rho_{L2}(r) R_{L2}(t,p) = 
!  = \sum_{L1,L3} w_{L1}^{*}(r) Y_{L1}^{*}(t,p) * w_{L3}(r) Y_{L3} (t,p)
! rho_{L2}(r) = \sum_{L1,L3} w_{L1}^{*}(r) w_{L3}(r)   <Y_{L1} |R_{L2}| Y_{L3}>
allocate(f1(s_nr_min,lmmaxwan,nspinor))
allocate(f2(s_nr_min,lmmaxwan,nspinor))
! rearrange wanlm in memory
do ispn=1,nspinor
  do lm=1,lmmaxwan
    f1(:,lm,ispn)=wanlm(lm,1:s_nr_min,ispn)
  enddo
enddo
f2=0.d0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(lm1,lm2,lm3,zt1,ispn)
do lmloc=1,lmmaxwanloc
  lm2=mpi_grid_map2(lmmaxwan,dims=(/dim_k,dim2/),loc=lmloc)
  do lm1=1,lmmaxwan
    do lm3=1,lmmaxwan
      zt1=gauntyry(lm2l(lm1),lm2l(lm2),lm2l(lm3),&
                   lm2m(lm1),lm2m(lm2),lm2m(lm3))
      if (abs(zt1).gt.1d-12) then
        do ispn=1,nspinor
          f2(:,lm2,ispn)=f2(:,lm2,ispn)+&
            dreal(dconjg(f1(:,lm1,ispn))*f1(:,lm3,ispn)*zt1)
        enddo !ispn
      endif
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
call mpi_grid_reduce(f2(1,1,1),s_nr_min*lmmaxwan*nspinor,all=.true.)
do ispn=1,nspinor
  do lm=1,lmmaxwan
    rholm(lm,:,ispn)=f2(:,lm,ispn)
  enddo
  totrholm(:,:)=totrholm(:,:)+rholm(:,:,ispn)
enddo
deallocate(f1,f2)
! norm of total charge density
t1=0.d0
do ir=1,s_nr_min
  t1=t1+totrholm(1,ir)*s_rw(ir)
enddo
wanprop(wp_normrho)=t1*fourpi*y00
! estimate the quadratic spread <r^2>-<r>^2
x2=0.d0
x=0.d0
! TODO: check this formulas
! Ry=-\frac{1}{2} \sqrt{\frac{3}{\pi }} \sin (t) \sin (p)
! Rz= \frac{1}{2} \sqrt{\frac{3}{\pi }} \cos (t)
! Rx=-\frac{1}{2} \sqrt{\frac{3}{\pi }} \sin (t) \cos (p) 
do ir=1,s_nr_min
  x(1)=x(1)-2.d0*s_r(ir)*sqrt(pi/3.d0)*totrholm(4,ir)*s_rw(ir)
  x(2)=x(2)-2.d0*s_r(ir)*sqrt(pi/3.d0)*totrholm(2,ir)*s_rw(ir)
  x(3)=x(3)+2.d0*s_r(ir)*sqrt(pi/3.d0)*totrholm(3,ir)*s_rw(ir)
  x2=x2+2.d0*(s_r(ir)**2)*sqrt(pi)*totrholm(1,ir)*s_rw(ir)
enddo
wanprop(wp_spread)=x2-dot_product(x,x)
wanprop(wp_spread_x)=x(1)
wanprop(wp_spread_y)=x(2)
wanprop(wp_spread_z)=x(3)
! compute Hartree potential
vhalm=0.d0
do lmloc=1,lmmaxwanloc
  lm=mpi_grid_map2(lmmaxwan,dims=(/dim_k,dim2/),loc=lmloc)
  l=lm2l(lm)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(jr,t1)
  do ir=1,s_nr_min
    t1=0.d0
    do jr=1,ir
      t1=t1+(s_r(jr)**l/s_r(ir)**(l+1))*totrholm(lm,jr)*s_rw(jr)
    enddo
    do jr=ir+1,s_nr_min
      t1=t1+(s_r(ir)**l/s_r(jr)**(l+1))*totrholm(lm,jr)*s_rw(jr)
    enddo
    vhalm(lm,ir)=t1*fourpi/(2*l+1)
  enddo
!$OMP END PARALLEL DO
enddo
call mpi_grid_reduce(vhalm(1,1),lmmaxwan*s_nr_min,all=.true.)
! compute XC potential
if (spinpol) then
  call xcifc(xctype,n=s_ntp*s_nr_min,rhoup=rhotp(1,1,1),rhodn=rhotp(1,1,2),&
    ex=extp,ec=ectp,vxup=vxtp(1,1,1),vxdn=vxtp(1,1,2),vcup=vctp(1,1,1),&
    vcdn=vctp(1,1,2))
else
  call xcifc(xctype,n=s_ntp*s_nr_min,rho=rhotp,ex=extp,ec=ectp,vx=vxtp,vc=vctp)
endif
! save exchange and correlation energies
!if (.true.) then
!  call sic_rbsht(s_nr_min,extp,exclm) 
!  wanprop(wp_ex)=0.d0
!  do ir=1,s_nr_min
!    wanprop(wp_ex)=wanprop(wp_ex)+&
!      ddot(lmmaxwan,totrholm(1,ir),1,exclm(1,ir),1)*s_rw(ir)
!  enddo
!  call sic_rbsht(s_nr_min,ectp,exclm) 
!  wanprop(wp_ec)=0.d0
!  do ir=1,s_nr_min
!    wanprop(wp_ec)=wanprop(wp_ec)+&
!      ddot(lmmaxwan,totrholm(1,ir),1,exclm(1,ir),1)*s_rw(ir)
!  enddo
!endif
! save XC energy density in extp
if (sicec) extp(:,:)=extp(:,:)+ectp(:,:)
! expand in real spherical harmonics
call sic_rbsht(s_nr_min,extp,exclm) 
! save XC potential in vxtp and expand in real spherical harmonics   
do ispn=1,nspinor
  if (sicvc) vxtp(:,:,ispn)=vxtp(:,:,ispn)+vctp(:,:,ispn)
  call sic_rbsht(s_nr_min,vxtp(1,1,ispn),vxclm(1,1,ispn))
enddo
! compute vha=<V_h|rho>
wanprop(wp_vha)=0.d0
do ir=1,s_nr_min
  wanprop(wp_vha)=wanprop(wp_vha)+&
    ddot(lmmaxwan,totrholm(1,ir),1,vhalm(1,ir),1)*s_rw(ir)
enddo
! compute exc=<E_xc|rho>
wanprop(wp_exc)=0.d0
do ir=1,s_nr_min
  wanprop(wp_exc)=wanprop(wp_exc)+&
    ddot(lmmaxwan,totrholm(1,ir),1,exclm(1,ir),1)*s_rw(ir)
enddo
! compute vxc=<W_n|V_xc|W_n>; in the collinear case this is 
!  \sum_{\sigma} <V_xc^{\sigma}|rho${\sigma}>
wanprop(wp_vxc)=0.d0
do ispn=1,nspinor
  do ir=1,s_nr_min
    wanprop(wp_vxc)=wanprop(wp_vxc)+&
      ddot(lmmaxwan,rholm(1,ir,ispn),1,vxclm(1,ir,ispn),1)*s_rw(ir)
  enddo
enddo
! compute <V_n|rho>
wanprop(wp_vsic)=wanprop(wp_vha)+wanprop(wp_vxc)
! add Hartree potential to XC
do ispn=1,nspinor
  vxclm(:,:,ispn)=vxclm(:,:,ispn)+vhalm(:,:)
enddo
!if (.true.) then
!  call sic_write_pot(n,vxclm)
!endif
deallocate(rhotp,rholm,totrholm)
deallocate(vhalm,extp,ectp,exclm,vxtp,vctp)
! multiply Wannier function with potential and change sign
call timer_start(t_sic_wvprod)
allocate(f1(s_nr_min,lmmaxwan,nspinor))
allocate(f2(s_nr_min,lmmaxwan,nspinor))
allocate(f3(s_nr_min,lmmaxwan,nspinor))
! rearrange arrays in memory
do ispn=1,nspinor
  do lm=1,lmmaxwan
    f1(:,lm,ispn)=wanlm(lm,1:s_nr_min,ispn)
    f2(:,lm,ispn)=vxclm(lm,:,ispn)
  enddo
enddo
f3=zzero
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(lm1,lm2,lm3,zt1,ispn)
do lmloc=1,lmmaxwanloc
  lm3=mpi_grid_map2(lmmaxwan,dims=(/dim_k,dim2/),loc=lmloc)
  do lm1=1,lmmaxwan
    do lm2=1,lmmaxwan
      zt1=gauntyry(lm2l(lm1),lm2l(lm2),lm2l(lm3),&
               lm2m(lm1),lm2m(lm2),lm2m(lm3))
      if (abs(zt1).gt.1d-12) then
        do ispn=1,nspinor
          f3(:,lm3,ispn)=f3(:,lm3,ispn)-f1(:,lm1,ispn)*f2(:,lm2,ispn)*zt1
        enddo !ispn
      endif
    enddo
  enddo
enddo
!$OMP END PARALLEL DO
call mpi_grid_reduce(f3(1,1,1),s_nr_min*lmmaxwan*nspinor,all=.true.) 
wvlm=zzero
do ispn=1,nspinor
  do lm=1,lmmaxwan
    wvlm(lm,1:s_nr_min,ispn)=f3(1:s_nr_min,lm,ispn)
  enddo
enddo
deallocate(f1,f2,f3,vxclm)
call timer_stop(t_sic_wvprod)
end subroutine   

subroutine sic_write_pot(n,vxclm)
use modmain
use mod_sic
use mod_hdf5
implicit none
integer, intent(in) :: n
real(8), intent(in) :: vxclm(lmmaxwan,s_nr_min,nspinor)
!
character*100 fname
!
write(fname,'("vlm__",I4.4,".hdf5")')n
call hdf5_create_file(trim(fname))  
call hdf5_write(trim(fname),"/","lmmaxwan",lmmaxwan)
call hdf5_write(trim(fname),"/","s_nr_min",s_nr_min)
call hdf5_write(trim(fname),"/","nspinor",nspinor)
call hdf5_write(trim(fname),"/","s_nr",s_nr)
call hdf5_write(trim(fname),"/","vlm",vxclm(1,1,1),&
  (/lmmaxwan,s_nr_min,nspinor/))
call hdf5_write(trim(fname),"/","s_r",s_r(1),(/s_nr/))
return
end subroutine







