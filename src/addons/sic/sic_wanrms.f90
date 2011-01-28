subroutine sic_wanrms(n,wantp,wanlm,wanrms)
use modmain
use mod_sic
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(in) :: wantp(s_ntp,s_nr,nspinor)
complex(8), intent(in) :: wanlm(lmmaxwan,s_nr,nspinor)
real(8), intent(out) :: wanrms(3)
! local variables
integer nr,ntp,ir,itp,ispn
real(8) t1
complex(8) zt1,wanval(nspinor)
real(8), allocatable :: tp(:,:),vrc(:,:,:)
complex(8), allocatable :: rhotp(:,:)
complex(8), allocatable :: rholm(:,:)
complex(8), allocatable :: wantpc(:,:,:)
complex(8), allocatable :: rhotpc(:,:)
!
nr=300
ntp=200

allocate(rhotp(s_ntp,s_nr))
allocate(rholm(lmmaxwan,s_nr))
allocate(wantpc(ntp,nr,nspinor))
allocate(rhotpc(ntp,nr))
allocate(tp(2,ntp))
allocate(vrc(3,ntp,nr))
call sphcover(ntp,tp)
! resolve Wannier function on a coarse mesh
do ir=1,nr
  do itp=1,ntp
    vrc(:,itp,ir)=(/sin(tp(1,itp))*cos(tp(2,itp)),&
                    sin(tp(1,itp))*sin(tp(2,itp)),&
                    cos(tp(1,itp))/)*ir*sic_wan_cutoff/nr/2.d0
    call s_get_wanval(n,vrc(1,itp,ir),wanval)
    wantpc(itp,ir,:)=wanval(:)
  enddo
enddo
call mpi_grid_reduce(wantpc(1,1,1),ntp*nr*nspinor,dims=(/dim_k/),all=.true.)
! compute RMS for Wannier functions
t1=0.d0
do ir=1,nr
  do itp=1,ntp
    do ispn=1,nspinor
      zt1=wantpc(itp,ir,ispn)-s_func_val(vrc(1,itp,ir),wanlm(1,1,ispn))
      t1=t1+abs(zt1)**2
    enddo
  enddo
enddo
wanrms(1)=sqrt(t1/nr/ntp)
! compute charge density
rhotp=zzero
do ispn=1,nspinor
  rhotp(:,:)=rhotp(:,:)+zone*abs(wantp(:,:,ispn))**2
enddo
! convert to spherical harmonics
call zgemm('T','N',lmmaxwan,s_nr,s_ntp,zone,s_ylmb,s_ntp,rhotp,&
  s_ntp,zzero,rholm,lmmaxwan)
! compute charge density on a coarse mesh
rhotpc=zzero
do ispn=1,nspinor
  rhotpc(:,:)=rhotpc(:,:)+zone*abs(wantpc(:,:,ispn))**2
enddo
! compute RMS for charge density
t1=0.d0
do ir=1,nr
  do itp=1,ntp
    zt1=rhotpc(itp,ir)-s_func_val(vrc(1,itp,ir),rholm)
    t1=t1+abs(zt1)**2
  enddo
enddo
wanrms(2)=sqrt(t1/nr/ntp)
! compute rho^{1/3}
rhotp(:,:)=rhotp(:,:)**(1/3.d0)
! convert to spherical harmonics
call zgemm('T','N',lmmaxwan,s_nr,s_ntp,zone,s_ylmb,s_ntp,rhotp,&
  s_ntp,zzero,rholm,lmmaxwan)
rhotpc(:,:)=rhotpc(:,:)**(1/3.d0)
! compute RMS for rho^{1/3}
t1=0.d0
do ir=1,nr
  do itp=1,ntp
    zt1=rhotpc(itp,ir)-s_func_val(vrc(1,itp,ir),rholm)
    t1=t1+abs(zt1)**2
  enddo
enddo
wanrms(3)=sqrt(t1/nr/ntp)
deallocate(rhotp,rholm,wantpc,rhotpc,tp,vrc)
return
end
