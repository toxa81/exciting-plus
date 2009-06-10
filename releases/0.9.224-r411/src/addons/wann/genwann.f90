subroutine genwann(ik,evecfv,evecsv)
use modmain
implicit none
! arguments
integer, intent(in) :: ik
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(in) :: evecsv(nstsv,nstsv)
! local variables
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfsvmt(:,:,:,:,:)
complex(8), allocatable :: wfsvit(:,:,:)
complex(8), allocatable :: wann_unkmt_new(:,:,:,:,:)
complex(8), allocatable :: wann_unkit_new(:,:,:)
integer j,n
integer, external :: ikglob
! allocate arrays
allocate(wfsvmt(lmmaxvr,nrfmax,natmtot,nspinor,nstsv))
allocate(wfsvit(ngkmax,nspinor,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
call match(ngk(1,ikglob(ik)),gkc(1,1,ik),tpgkc(1,1,1,ik),sfacgk(1,1,1,ik),apwalm)
! generate second-varioational wave-functions
call genwfsvmt(lmaxvr,lmmaxvr,ngk(1,ikglob(ik)),evecfv,evecsv,apwalm,wfsvmt)
call genwfsvit(ngk(1,ikglob(ik)),evecfv,evecsv,wfsvit)
! calculate WF expansion coefficients
call genwann_c(ikglob(ik),evalsv(1,ikglob(ik)),wfsvmt,wann_c(1,1,ik))
! compute Bloch-sums of Wannier functions
allocate(wann_unkmt_new(lmmaxvr,nrfmax,natmtot,nspinor,nwann))
allocate(wann_unkit_new(ngkmax,nspinor,nwann))
wann_unkmt_new=zzero
wann_unkit_new=zzero
do n=1,nwann
  do j=1,nstsv
    wann_unkmt_new(:,:,:,:,n)=wann_unkmt_new(:,:,:,:,n) + &
      wfsvmt(:,:,:,:,j)*wann_c(n,j,ik)
    wann_unkit_new(:,:,n)=wann_unkit_new(:,:,n) + &
      wfsvit(:,:,j)*wann_c(n,j,ik)
  enddo
enddo
wann_unkmt(:,:,:,:,:,ik)=wann_unkmt_new(:,:,:,:,:)
wann_unkit(:,:,:,ik)=wann_unkit_new(:,:,:)
deallocate(wfsvmt,wfsvit,apwalm,wann_unkmt_new,wann_unkit_new)
return
end


