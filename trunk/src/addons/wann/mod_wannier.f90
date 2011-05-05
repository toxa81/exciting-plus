module mod_wannier
implicit none

!------------------!
!     Wannier      !
!------------------!
logical wannier
integer wann_natom
integer wann_norbgrp
integer wann_ntype
logical wann_add_poco
integer, allocatable :: wann_norb(:)
integer, allocatable :: wann_iorb(:,:,:)
integer, allocatable :: wann_iprj(:,:)
real(8), allocatable :: wann_eint(:,:)
real(8), allocatable :: wann_v(:)

integer nwantot
integer, allocatable :: wan_info(:,:)
integer, allocatable :: nwannias(:)
! positions of WFs
real(8), allocatable :: wanpos(:,:)
! magnetic moment of WFs
real(8), allocatable :: wanmom(:)
! expansion coefficients of Wannier functions over spinor Bloch eigen-functions  
complex(8), allocatable :: wann_c(:,:,:)
! Bloch-sums of WF
complex(8), allocatable :: wann_unkmt(:,:,:,:,:,:)
complex(8), allocatable :: wann_unkit(:,:,:,:)

! H(k) in WF basis
complex(8), allocatable :: wann_h(:,:,:)
! e(k) of WF H(k) (required for band-sctructure plot only)
real(8), allocatable :: wann_e(:,:)
! momentum operator in WF basis
complex(8), allocatable :: wann_p(:,:,:,:)

real(8), allocatable :: wann_ene(:)
real(8), allocatable :: wann_occ(:)

real(8) zero3d(3)
real(8) bound3d(3,3)
integer nrxyz(3)
integer nwfplot
integer firstwf
logical wannier_lc
! linear combination of projected orbitals for Wannier functions
data wannier_lc/.false./
! number of linear combinations; this defines the new number of Wannier functions
integer nwanlc
! number of atom-centerd orbitals to make each mulipole orbital
integer, allocatable :: wanlc_norb(:)
! index and tranlation of each atom-cented orbital
integer, allocatable :: wanlc_iorb(:,:,:)
! weight of each atom-cented orbital
real(8), allocatable :: wanlc_iorb_alpha(:,:)
! maximum number of orbials to make a mulipole orbital
integer, parameter :: wanlcmax=100 

integer nwann_h
integer, allocatable :: iwann_h(:)

logical wannier_soft_eint
data wannier_soft_eint/.false./
real(8), allocatable :: wannier_soft_eint_w1(:)
real(8), allocatable :: wannier_soft_eint_w2(:)
real(8), allocatable :: wannier_soft_eint_e1(:)
real(8), allocatable :: wannier_soft_eint_e2(:)
real(8) wannier_min_prjao
data wannier_min_prjao/-0.1d0/

logical ldisentangle
data ldisentangle/.false./

logical lrespwffilter
data lrespwffilter/.false./
integer nwfch
integer wfch(2,200)

! wannier_prjao=0 (default) : project to local orbital
! wannier_prjao=1 : project to f(x)=(1+cos(Pi*x/R))
integer wannier_prjao
data wannier_prjao/0/

integer, allocatable :: wann_err_k(:)

integer, allocatable :: wannier_prjlo(:,:)



type wannier_transitions
! number of taken Wannier functions
  integer :: nwan
! global index of taken Wannier functions
  integer, allocatable :: iwan(:)
! mapping from global to local index
  integer, allocatable :: idxiwan(:) 
! total number of Wannier transitions (total number of <m| |nT> bra-kets)
  integer :: nwt
! i-th Wannier transition
  integer, allocatable :: iwt(:,:)
! mapping from {m,n,T} to global index
  integer, allocatable :: iwtidx(:,:,:,:,:) 
! translation limits
  integer :: tlim(2,3)  
! minimal distance
  real(8) :: mindist
! maximum distance
  real(8) :: maxdist
! number of all encountered translations
  integer :: ntr
! list of translations
  integer, allocatable :: vtr(:,:)
! flag that indicates if transition between m,n Wannier functions is in list
  integer, allocatable :: wt(:,:)
end type wannier_transitions

contains

subroutine genwantran(twantran,mindist,maxdist,waninc,allwt,diagwt)
use mod_addons
implicit none
! arguments
type(wannier_transitions), intent(out) :: twantran
real(8), intent(in) :: mindist
real(8), intent(in) :: maxdist
integer, optional, intent(in) :: waninc(nwantot)
logical, optional, intent(in) :: allwt
logical, optional, intent(in) :: diagwt
! local variables
integer m,n,j,nwtmax,i,k,t(3),ias,jas
logical ladd,allwt_,diagwt_
integer nwan,nwt,ntr,tlim(2,3)
integer, allocatable :: iwan(:)
integer, allocatable :: idxiwan(:)
integer, allocatable :: iwt(:,:)
integer, allocatable :: wt(:,:)
integer, allocatable :: vtr(:,:)
logical, external :: wann_diel
real(8), parameter :: epswfocc=1d-8

allwt_=.false.
if (present(allwt)) allwt_=allwt
diagwt_=.false.
if (present(diagwt)) diagwt_=diagwt

! make list of included Wannier functions
allocate(iwan(nwantot))
if (present(waninc)) then
  nwan=0
  do n=1,nwantot
    if (waninc(n).ne.0) then
      nwan=nwan+1
      iwan(nwan)=n
    endif   
  enddo
else
  nwan=nwantot
  do n=1,nwantot
    iwan(n)=n
  enddo
endif
twantran%nwan=nwan
allocate(twantran%iwan(nwan))
twantran%iwan(1:nwan)=iwan(1:nwan)
! mapping from global to local index
allocate(idxiwan(nwantot))
idxiwan=-1
do i=1,nwan
  n=iwan(i)
  idxiwan(n)=i
enddo
allocate(twantran%idxiwan(nwantot))
twantran%idxiwan=idxiwan
! get nearest neighbours
twantran%mindist=mindist
twantran%maxdist=maxdist
call getnghbr(mindist,maxdist)
! get maximum possible number of WF transitions
nwtmax=0
do n=1,nwantot
  ias=wan_info(1,n)
  do i=1,nnghbr(ias)
    jas=inghbr(1,i,ias)
    nwtmax=nwtmax+nwannias(jas)
  enddo
enddo
! make list of Wannier transitions
allocate(iwt(5,nwtmax))
allocate(wt(nwantot,nwantot))
wt=0
nwt=0
do i=1,nwan
  m=iwan(i)
  ias=wan_info(1,m)
  do k=1,nnghbr(ias)
    do j=1,nwan
      n=iwan(j)
      jas=wan_info(1,n)
      if (jas.eq.inghbr(1,k,ias)) then
        ladd=.false.
        if (diagwt_) then
          if (m.eq.n) ladd=.true.  
        else
! for integer occupancy numbers take only transitions between occupied and empty bands
          if (wann_diel().and.(abs(wann_occ(m)-wann_occ(n)).gt.epswfocc)) ladd=.true.
! for fractional occupancies or other cases take all transitions
          if (.not.wann_diel().or.allwt_) ladd=.true.
        endif
        if (ladd) then
          nwt=nwt+1
          iwt(1,nwt)=m
          iwt(2,nwt)=n
          iwt(3:5,nwt)=inghbr(3:5,k,ias)
          wt(m,n)=1
        endif
      endif
    enddo !j1
  enddo !i
enddo !j
twantran%nwt=nwt
allocate(twantran%iwt(5,nwt))
twantran%iwt(:,1:nwt)=iwt(:,1:nwt)
allocate(twantran%wt(nwantot,nwantot))
twantran%wt(:,:)=wt(:,:)
! get translation limits
tlim=0
if (nwt.gt.0) then
  do i=1,3
    tlim(1,i)=minval(iwt(2+i,1:nwt))
    tlim(2,i)=maxval(iwt(2+i,1:nwt))
  enddo
endif
twantran%tlim(:,:)=tlim(:,:)
! generate {m,n,t} -> global index mapping
allocate(twantran%iwtidx(nwantot,nwantot,tlim(1,1):tlim(2,1),&
  tlim(1,2):tlim(2,2),tlim(1,3):tlim(2,3)))
twantran%iwtidx=-1
do i=1,nwt
  m=iwt(1,i)
  n=iwt(2,i)
  t=iwt(3:5,i)
  twantran%iwtidx(m,n,t(1),t(2),t(3))=i
enddo
! get list of encountered translations
allocate(vtr(3,nwt))
ntr=0
do i=1,nwt
  t(:)=iwt(3:5,i)
  ladd=.true.
  do j=1,ntr
    if (all(vtr(:,j).eq.t(:))) ladd=.false.
  enddo
  if (ladd) then
    ntr=ntr+1
    vtr(:,ntr)=t(:)
  endif
enddo
twantran%ntr=ntr
allocate(twantran%vtr(3,ntr))
twantran%vtr(:,1:ntr)=vtr(:,1:ntr)
deallocate(iwan)
deallocate(idxiwan)
deallocate(iwt)
deallocate(wt)
deallocate(vtr)
return
end subroutine

subroutine deletewantran(twantran)
implicit none
type(wannier_transitions), intent(inout) :: twantran

if (allocated(twantran%iwan)) deallocate(twantran%iwan)
if (allocated(twantran%idxiwan)) deallocate(twantran%idxiwan)
if (allocated(twantran%iwt)) deallocate(twantran%iwt)
if (allocated(twantran%iwtidx)) deallocate(twantran%iwtidx)
if (allocated(twantran%vtr)) deallocate(twantran%vtr)
if (allocated(twantran%wt)) deallocate(twantran%wt)
return
end subroutine

subroutine printwantran(twantran)
implicit none
type(wannier_transitions), intent(in) :: twantran
integer i

write(*,*)'twantran%nwt=',twantran%nwt
do i=1,twantran%nwt
  write(*,*)twantran%iwt(:,i)
enddo
return
end subroutine

end module
