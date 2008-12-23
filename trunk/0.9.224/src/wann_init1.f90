subroutine wann_init1
use modmain
implicit none

real(8) r(3),r1(3),t(3),d
real(8) bound3d(3,3),orig3d(3),zero3d(3)
real(8) bound2d(3,2),orig2d(3)
complex(4), allocatable :: wf(:,:)

integer ntr(3),i,ivec,nrxyz(3),nrtot
integer i1,i2,i3,ir,n,ierr
integer ik,ispn,istfv,j,ig
complex(8), allocatable :: evecfv(:,:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: evec(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: bcoeff(:,:,:,:,:,:,:)
complex(8), allocatable :: ccoeff(:,:,:)
complex(8), allocatable :: wfsvitloc(:,:,:)
complex(8), allocatable :: wfsvmtloc(:,:,:,:,:)
complex(8) zt1
complex(8), allocatable :: zt2(:,:,:,:,:)
complex(8), allocatable :: zt3(:,:,:,:)
character*8 fname
real(8) x(2),alph
logical, parameter :: wf3d=.true.
integer tlim(2,3)
integer nwfplot,firstwf
integer, external :: ikglob

! read density and potentials from file
call readstate
! read Fermi energy from file
call readfermi
! generate the core wavefunctions and densities
call gencore
! find the new linearisation energies
call linengy
! generate the APW radial functions
call genapwfr
! generate the local-orbital radial functions
call genlofr

call geturf

nwfplot=5
firstwf=1

allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(evecsv(nstsv,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wfsvmtloc(lmmaxvr,nrfmax,natmtot,nstsv,nspinor))
allocate(wfsvitloc(nmatmax,nstsv,nspinor))

allocate(wann_unkmt(lmmaxvr,nrfmax,natmtot,wf_dim,nspinor,nkpt))
allocate(wann_unkit(nmatmax,wf_dim,nspinor,nkpt))
wann_unkmt=dcmplx(0.d0,0.d0)
wann_unkit=dcmplx(0.d0,0.d0)

! read and transform eigen-vectors
do ik=1,nkpt
  call getevecfv(vkl(1,ik),vgkl(1,1,1,ik),evecfv)
  call getevecsv(vkl(1,ik),evecsv)
  call getwfc(ik,wfc(1,1,1,ik))
  call match(ngk(1,ik),gkc(1,1,ik),tpgkc(1,1,1,ik),sfacgk(1,1,1,ik),apwalm)
  call genwfsvmt(lmaxvr,lmmaxvr,ngk(1,ik),evecfv,evecsv,apwalm,wfsvmtloc)
  call genwfsvit(ngk(1,ik),evecfv,evecsv,wfsvitloc)
  do n=1,wf_dim
    do i=1,nstfv
      wann_unkmt(:,:,:,n,1,ik)=wann_unkmt(:,:,:,n,1,ik)+wfsvmtloc(:,:,:,i,1)*wfc(n,i,1,ik)
      wann_unkit(:,n,1,ik)=wann_unkit(:,n,1,ik)+wfsvitloc(:,i,1)*wfc(n,i,1,ik)
    enddo
  enddo
enddo !ik

return
end

subroutine wann_unk(n,vpl,vrc,val)
use modmain
implicit none
integer, intent(in) :: n
real(8), intent(in) :: vpl(3)
real(8), intent(in) :: vrc(3)
real(8), intent(out) :: val(2)



return
end


