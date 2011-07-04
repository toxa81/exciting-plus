subroutine sic_gensmesh
use modmain
use mod_sic
implicit none
integer itp,lm
real(8) a,tp(2)
integer nt,np,it,ip,l,m,lm1
real(8) t1
complex(8) zt1
complex(8), allocatable :: clm(:)
complex(8), allocatable :: zftp(:)
!
!call gensmesh_icos
call gensmesh_cubed

if (allocated(s_rlmf)) deallocate(s_rlmf)
allocate(s_rlmf(lmmaxwan,s_ntp))
if (allocated(s_ylmf)) deallocate(s_ylmf)
allocate(s_ylmf(lmmaxwan,s_ntp))
if (allocated(s_rlmb)) deallocate(s_rlmb)
allocate(s_rlmb(s_ntp,lmmaxwan))
if (allocated(s_ylmb)) deallocate(s_ylmb)
allocate(s_ylmb(s_ntp,lmmaxwan))
! generate spherical harmonics
do itp=1,s_ntp
  call genrlm(lmaxwan,s_tp(1,itp),s_rlmf(1,itp))
  call genylm(lmaxwan,s_tp(1,itp),s_ylmf(1,itp))
  do lm=1,lmmaxwan
    s_rlmb(itp,lm)=s_rlmf(lm,itp)*s_tpw(itp)
    s_ylmb(itp,lm)=dconjg(s_ylmf(lm,itp))*s_tpw(itp) 
  enddo
enddo

call test_smesh

!call gensmesh_lebedev
!call gensmesh_uniform
!call gensmesh_cubed
!call gensmesh_symcrys
!call gensmesh_healpix
!call gensmesh_icos

!call bstop


!call bstop
!call nfsft_init(lmaxwan,s_ntp,stp)
return
end subroutine

subroutine test_smesh
use modmain
use mod_sic
implicit none
integer i,lm,lm1,itp
real(8) t1,t2
complex(8), allocatable :: zflm(:),zftp(:),zflm_orig(:)
! test 1: reconstruction of random vector
allocate(zflm_orig(lmmaxwan))
allocate(zflm(lmmaxwan))
allocate(zftp(s_ntp))
do i=0,2
  zflm_orig=zzero
  do lm=1,lmmaxwan
    call random_number(t1); t1=t1-0.5d0
    call random_number(t2); t2=t2-0.5d0
    zflm_orig(lm)=dcmplx(t1,t2)*(10**i)
  enddo
  zftp=zzero
  do itp=1,s_ntp
    do lm=1,lmmaxwan
      zftp(itp)=zftp(itp)+zflm_orig(lm)*s_ylmf(lm,itp)
    enddo
  enddo
  call sic_zbsht(1,zftp,zflm)
  t1=0.d0
  do lm=1,lmmaxwan
    t1=t1+abs(zflm(lm)-zflm_orig(lm))
  enddo
  if (mpi_grid_root()) then
    write(*,'("[test_smesh] complex random vector reconstruction error : ",G18.10)')t1
  endif
enddo
! test 2: backward transformation of spherical harmonics
if (sic_debug_level.ge.3) then 
  do lm=1,lmmaxwan
    zftp(:)=s_ylmf(lm,:) 
    call sic_zbsht(1,zftp,zflm)
    t1=0.d0
    do lm1=1,lmmaxwan
      if (lm1.eq.lm) zflm(lm1)=zflm(lm1)-zone
      t1=t1+abs(zflm(lm1))
    enddo
    t2=0.d0
    zflm=zzero
    do lm1=1,lmmaxwan
      do itp=1,s_ntp
        zflm(lm1)=zflm(lm1)+s_ylmb(itp,lm1)*zftp(itp)
      enddo
      if (lm1.eq.lm) zflm(lm1)=zflm(lm1)-zone
      t2=t2+abs(zflm(lm1))
    enddo
    write(*,'("lm=",I4,"    sic_zbsht error : ",G18.10,&
      &"  single summation error : ",G18.10)')lm,t1,t2
  enddo
endif
deallocate(zflm_orig)
deallocate(zflm)
deallocate(zftp)
! test 3: check weights
!do l=0,lmaxwan
!  t1=0.d0
!  do m=-l,l
!    lm=idxlm(l,m)
!    zt1=zzero
!    do itp=1,s_ntp
!      zt1=zt1+s_ylmb(itp,lm)
!    enddo
!    t1=t1+abs(zt1)**2
!  enddo
!  write(500,*)l,t1
!enddo
return
end subroutine








