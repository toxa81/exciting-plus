
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: seceqnfv
! !INTERFACE:
subroutine seceqnfv(ik,nmatp,ngp,igpig,vgpc,apwalm,evalfv,evecfv)
! !USES:
use modmain
use mod_seceqn
! !INPUT/OUTPUT PARAMETERS:
!   nmatp  : order of overlap and Hamiltonian matrices (in,integer)
!   ngp    : number of G+k-vectors for augmented plane waves (in,integer)
!   igpig  : index from G+k-vectors to G-vectors (in,integer(ngkmax))
!   vgpc   : G+k-vectors in Cartesian coordinates (in,real(3,ngkmax))
!   apwalm : APW matching coefficients
!            (in,complex(ngkmax,apwordmax,lmmaxapw,natmtot))
!   evalfv : first-variational eigenvalues (out,real(nstfv))
!   evecfv : first-variational eigenvectors (out,complex(nmatmax,nstfv))
! !DESCRIPTION:
!   Solves the secular equation,
!   $$ (H-\epsilon O)b=0, $$
!   for the all the first-variational states of the input $k$-point.
!
! !REVISION HISTORY:
!   Created March 2004 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ik
integer, intent(in) :: nmatp
integer, intent(in) :: ngp
integer, intent(in) :: igpig(ngkmax)
real(8), intent(in) :: vgpc(3,ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
real(8), intent(out) :: evalfv(nstfv)
complex(8), intent(out) :: evecfv(nmatmax,nstfv)
! local variables
integer is,ia,i,m,np,info,nb,lwork
real(8) v(1),vl,vu
real(8) ts0,ts1
! allocatable arrays
integer, allocatable :: iwork(:)
integer, allocatable :: ifail(:)
real(8), allocatable :: w(:)
real(8), allocatable :: rwork(:)
complex(8), allocatable :: h(:),h1(:)
complex(8), allocatable :: o(:),o1(:)
complex(8), allocatable :: work(:)
complex(8) zt1
complex(4) ct1
character*100 fname
logical, parameter :: packed=.false.
integer, external :: ilaenv
!call lapw_seceqnfv(ngp,nmatp,nmatp,igpig,vgpc,veffig,cfunig,apwalm,apwfr,&
!    apwdfr,haa,hloa,hlolo,oalo,ololo,h1,o1,evalfv,evecfv)
!return
if (packed) then
  np=(nmatp*(nmatp+1))/2
else
  np=nmatp*nmatp
endif  
!-----------------------------------------------!
!     Hamiltonian and overlap matrix set up     !
!-----------------------------------------------!
call timesec(ts0)
call timer_start(t_seceqnfv)
call timer_start(t_seceqnfv_setup)
allocate(h(np),o(np))
!allocate(h1(np),o1(np))
h(:)=zzero
o(:)=zzero
!h1(:)=zzero
!o1(:)=zzero
if (packed) then
! Hamiltonian
  do is=1,nspecies
    do ia=1,natoms(is)
      call hmlaa(.false.,is,ia,ngp,apwalm,v,h)
      call hmlalo(.false.,is,ia,ngp,apwalm,v,h)
      call hmllolo(.false.,is,ia,ngp,v,h)
    end do
  end do
  call hmlistl(.false.,ngp,igpig,vgpc,v,h)
! overlap
  do is=1,nspecies
    do ia=1,natoms(is)
      call olpaa(.false.,is,ia,ngp,apwalm,v,o)
      call olpalo(.false.,is,ia,ngp,apwalm,v,o)
      call olplolo(.false.,is,ia,ngp,v,o)
    end do
  end do
  call olpistl(.false.,ngp,igpig,v,o)
else
  call sethml(ngp,nmatp,vgpc,igpig,apwalm,h)
  call setovl(ngp,nmatp,igpig,apwalm,o)
  !call lapw_seceqnfv(ngp,nmatp,nmatp,igpig,vgpc,veffig,cfunig,apwalm,apwfr,&
  !  apwdfr,haa,hloa,hlolo,oalo,ololo,h1,o1,evalfv,evecfv)
  !write(*,*)sum(abs(h-h1))
  !write(*,*)sum(abs(o-o1))
endif
call timesec(ts1)
call timer_stop(t_seceqnfv_setup)
timemat=timemat+ts1-ts0
!------------------------------------!
!     solve the secular equation     !
!------------------------------------!
if (mpi_grid_root((/dim2/))) then
  call timesec(ts0)
  call timer_start(t_seceqnfv_diag)
  evecfv=zzero
  ! LAPACK 3.x call
  if (packed) then
    allocate(iwork(5*nmatp))
    allocate(ifail(nmatp))
    allocate(w(nmatp))
    allocate(rwork(7*nmatp))
    allocate(work(2*nmatp))
    call zhpgvx(1,'V','I','U',nmatp,h,o,vl,vu,1,nstfv,evaltol,m,w,evecfv,nmatmax, &
      work,rwork,iwork,ifail,info)
    evalfv(1:nstfv)=w(1:nstfv)
    deallocate(iwork,ifail,w,rwork,work)
    if (info.ne.0) then
      write(*,*)
      write(*,'("Error(seceqnfv): diagonalisation failed")')
      write(*,'(" ZHPGVX returned INFO = ",I8)') info
      if (info.gt.nmatp) then
        i=info-nmatp
        write(*,'(" The leading minor of the overlap matrix of order ",I8)') i
        write(*,'("  is not positive definite")')
        write(*,'(" Order of overlap matrix : ",I8)') nmatp
        write(*,*)
      end if
      call pstop
    end if
  else
    call diagzheg(nmatp,nstfv,nmatmax,evaltol,h,o,evalfv,evecfv)
  endif
  call timesec(ts1)
  call timer_stop(t_seceqnfv_diag)
endif
call mpi_grid_bcast(evecfv(1,1),nmatmax*nstfv,dims=(/dim2/))
call mpi_grid_bcast(evalfv(1),nstfv,dims=(/dim2/))
timefv=timefv+ts1-ts0
call timer_stop(t_seceqnfv)
deallocate(h,o)
!deallocate(h1,o1)
return
end subroutine
!EOC

