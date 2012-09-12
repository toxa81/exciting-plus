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

integer, allocatable :: wan_complex(:,:)

integer, parameter :: wi_atom = 1
integer, parameter :: wi_lm = 2
integer, parameter :: wi_spin = 3
integer, parameter :: wi_type = 4
integer, parameter :: wi_group = 5
integer, parameter :: wi_idxa = 6
integer, parameter :: wi_idxo = 7


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
! number of Wannier transitions with zero translation (on-site transitions)
  integer :: nwt0
! transition with zero translation (global index ranging from 1 to nwt)
  integer, allocatable :: iwt0(:)
! mapping from global index [1,nwt] to T=0 transitions [1,nwt0]
  integer, allocatable :: iwt0idx(:)
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

! for Madness way of resolving WFs: each thread of each MPI rank can compute 
!  value of the function at the arbitrary r-point 
integer m_ngvec
integer, allocatable :: m_ngknr(:)
integer, allocatable :: m_igkignr(:,:)
complex(8), allocatable :: m_wann_unkmt(:,:,:,:,:)
complex(8), allocatable :: m_wann_unkit(:,:,:)
real(8) m_wanpos(3)

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
integer nwan,nwt,ntr,tlim(2,3),nwt0
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
  &tlim(1,2):tlim(2,2),tlim(1,3):tlim(2,3)))
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
! get transitions with zero translation
nwt0=0
do i=1,twantran%nwt
  m=twantran%iwt(1,i)
  ias=wan_info(1,m)
  n=twantran%iwt(2,i)
  jas=wan_info(1,n)
  if (ias.eq.jas.and.all(twantran%iwt(3:5,i).eq.0)) nwt0=nwt0+1
enddo
twantran%nwt0=nwt0
allocate(twantran%iwt0(nwt0))
allocate(twantran%iwt0idx(nwt))
twantran%iwt0idx=-1
nwt0=0
do i=1,twantran%nwt
  m=twantran%iwt(1,i)
  ias=wan_info(1,m)
  n=twantran%iwt(2,i)
  jas=wan_info(1,n)
  if (ias.eq.jas.and.all(twantran%iwt(3:5,i).eq.0)) then
    nwt0=nwt0+1
    twantran%iwt0(nwt0)=i
    twantran%iwt0idx(i)=nwt0
  endif
enddo
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
if (allocated(twantran%iwt0)) deallocate(twantran%iwt0)
if (allocated(twantran%iwt0idx)) deallocate(twantran%iwt0idx)
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

! auxiliary subroutine to generate WF expansion coefficients over 
!  second-variational wave-functions 
subroutine wan_gencsv_aux(ikloc,evecfv,evecsv,evecfd)
use modmain
use mod_seceqn
implicit none
integer, intent(in) :: ikloc
complex(8), optional, intent(in) :: evecfv(nmatmax,nstfv)
complex(8), optional, intent(in) :: evecsv(nstsv,nstsv)
complex(8), optional, intent(in) :: evecfd(nspinor*nmatmax,nstsv)
!
integer ik,lmax,lmmax,n,j
real(8) t1
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evec(:,:)
complex(8), allocatable :: wfsvmt(:,:,:,:,:)
!
lmax=3
lmmax=(lmax+1)**2
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
allocate(apwalm(ngkmax,lmmaxapw,apwordmax,natmtot))
call genapwalm(ngk(1,ik),gkc(1,1,ikloc),tpgkc(1,1,1,ikloc),&
  &sfacgk(1,1,1,ikloc),apwalm)
allocate(wfsvmt(lmmax,nufrmax,natmtot,nspinor,nstsv))
if (present(evecfv).and.present(evecsv)) then
  allocate(evec(nspinor*nmatmax,nstsv))
  call evecsvfd(evecfv,evecsv,evec)
  call genwfsvc(lmax,lmmax,ngk(1,ik),nstsv,apwalm,evec,wfsvmt)
  deallocate(evec)
else
  call genwfsvc(lmax,lmmax,ngk(1,ik),nstsv,apwalm,evecfd,wfsvmt)
endif

call wan_gencsv(lmmax,vkc(1,ik),evalsv(1,ik),wfsvmt,wann_c(1,1,ikloc),&
  &wann_err_k(ikloc))

if (.false.) then
  write(*,'("WF expansion coefficients for k-point : ",I2)')ik
  do j=1,nstsv
    t1=0.d0
    do n=1,nwantot
      t1=t1+abs(wann_c(n,j,ikloc))
    enddo
    if (t1.gt.1d-10) then
      write(*,'(" j : ",I4,"  wann_c : ",255F8.4)')j,(abs(wann_c(n,j,ikloc)),n=1,nwantot)
    endif
  enddo
endif
deallocate(apwalm,wfsvmt)
return
end subroutine

subroutine wan_gencsv(lmmax,vpc,eval,wfsvmt,wanc,ierr)
use modmain
implicit none
integer, intent(in) :: lmmax
real(8), intent(in) :: vpc(3)
real(8), intent(in) :: eval(nstsv)
complex(8), intent(in) :: wfsvmt(lmmax,nufrmax,natmtot,nspinor,nstsv)
complex(8), intent(inout) :: wanc(nwantot,nstsv)
integer, optional, intent(out) :: ierr
!
integer n,j,ias,lm,ispn,itype,ilo,ierr_,prj_to_ylm
logical tbndint 
real(8) d1
logical, external :: bndint
!
wanc=zzero
do n=1,nwantot
  if (.not.wannier_lc) then
    ias=wan_info(1,n)
    lm=wan_info(2,n)
    ispn=wan_info(3,n)
    itype=wan_info(4,n)
    ilo=wannier_prjlo(wan_info(7,n),wan_info(6,n))
    prj_to_ylm=wan_complex(wan_info(7,n),wan_info(6,n))
    do j=1,nstsv
      tbndint=bndint(j,eval(j),wann_eint(1,itype),wann_eint(2,itype))
      if (tbndint) then
        call wan_genprjlo(ilo,prj_to_ylm,ias,lm,ispn,lmmax,wfsvmt(1,1,1,1,j),&
          wanc(n,j))
        if (wannier_soft_eint) then
          wanc(n,j)=wanc(n,j)*smoothstep(eval(j),wannier_soft_eint_e1(itype), &
            wannier_soft_eint_e2(itype),wannier_soft_eint_w1(itype),&
            wannier_soft_eint_w2(itype))
        endif
      endif
    enddo
  else
    stop "(wan_gencsv) linear combination of WFs has to be uncommented and tested"
  !  do i=1,wanlc_norb(n)
  !    d1=wanlc_iorb_alpha(i,n)
  !    iw=wanlc_iorb(1,i,n)
  !    itr(:)=wanlc_iorb(2:4,i,n)
  !    tr(:)=avec(:,1)*itr(1)+avec(:,2)*itr(2)+avec(:,3)*itr(3)
  !    ias=wan_info(1,iw)
  !    lm=wan_info(2,iw)
  !    ispn=wan_info(3,iw)
  !    itype=wan_info(4,iw)
  !    do j=1,nstsv
  !      if (bndint(j,e(j),wann_eint(1,itype),wann_eint(2,itype))) then
  !        call genprjao(ias,lm,ispn,j,wfsvmt,zt1)
  !        if (wannier_soft_eint) then
  !          zt1=zt1*orbwt(e(j),wannier_soft_eint_e1(itype),&
  !            wannier_soft_eint_e2(itype),wannier_soft_eint_w1(itype),&
  !            wannier_soft_eint_w2(itype))
  !        endif               
! !<psi_k(r)|g(r-T)>=<psi(r+T)|g(r)>=e^{-ikT}<psi(r)|g(r)>
  !        prjao(n,j)=prjao(n,j)+zt1*d1*exp(-zi*dot_product(vpc,tr))
  !      endif !bndint
  !    enddo
  !  enddo !i
  endif
enddo !n
! remove small contribution
do j=1,nstsv
  d1=0.d0
  do n=1,nwantot
    d1=d1+abs(wanc(n,j))**2
  enddo
  if (d1.lt.wannier_min_prjao) wanc(:,j)=zzero
enddo
!if (.false.) then
!  write(*,'("Total contribution of projected orbitals : ")')
!  do j=1,nstsv
!    d1=0.d0
!    do n=1,nwantot
!      d1=d1+abs(prjlo(n,j))**2
!    enddo
!    write(*,'("  band : ",I4,"  wt : ",F12.6)')j,d1
!  enddo
!endif
call wan_ort_k(wanc,ierr_)
!ierr_=0
if (present(ierr)) ierr=ierr_
return
end subroutine


subroutine wan_genprjlo(ilo,prj_to_ylm,ias,lm,ispn,lmmax,wfsvmt,prjlo)
use modmain
implicit none
integer, intent(in) :: ilo
integer, intent(in) :: prj_to_ylm
integer, intent(in) :: ias
integer, intent(in) :: lm
integer, intent(in) :: ispn
integer, intent(in) :: lmmax
complex(8), intent(in) :: wfsvmt(lmmax,nufrmax,natmtot,nspinor)
complex(8), intent(out) :: prjlo
!
integer l,m1,lm1,io1,ir,is,ic
real(8) fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)
!
! compute <psi_{ik}|phi_n>, where n={ias,lm,ispn} 
! |psi> is a spinor Bloch-function 
! |phi> is a valence local orbital
! 
l=lm2l(lm)
is=ias2is(ias)
ic=ias2ic(ias)
prjlo=zzero
if (prj_to_ylm.ne.0) then
  do io1=1,nufr(l,is)
    prjlo=prjlo+dconjg(wfsvmt(lm,io1,ias,ispn))*ufrp(l,io1,apword(l,is)+ilo,ic)
  enddo
else
  do m1=-l,l
    lm1=idxlm(l,m1)
    do io1=1,nufr(l,is)
  ! project to local orbital    
      if (wannier_prjao.eq.0) then
        prjlo=prjlo+dconjg(wfsvmt(lm1,io1,ias,ispn))*&
          &ufrp(l,io1,apword(l,is)+ilo,ic)*rylm_lps(lm,lm1,ias)
      endif
  ! project to f(x)=(1+cos(Pi*x/R))
      if (wannier_prjao.eq.1) then
        do ir=1,nrmt(is)
          fr(ir)=ufr(ir,l,io1,ic)*(1+cos(pi*spr(ir,is)/rmt(is)))*(spr(ir,is)**2)
        enddo
        call fderiv(-1,nrmt(is),spr(1,is),fr,gr,cf)
        prjlo=prjlo+dconjg(wfsvmt(lm1,io1,ias,ispn))*gr(nrmt(is))*rylm_lps(lm,lm1,ias)
      endif
    enddo !io1
  enddo !m
endif
return
end subroutine

subroutine wan_ort_k(wanc,ierr_)
use modmain
implicit none
complex(8), intent(inout) :: wanc(nwantot,nstsv)
integer, intent(out) :: ierr_
!
complex(8), allocatable :: s(:,:)
complex(8), allocatable :: wanc_ort(:,:)
integer m1,m2,j,ierr
integer, external :: hash
!
allocate(s(nwantot,nwantot))
allocate(wanc_ort(nwantot,nstsv))
! compute ovelap matrix
s=zzero
do m1=1,nwantot
  do m2=1,nwantot
    do j=1,nstsv
      s(m1,m2)=s(m1,m2)+dconjg(wanc(m1,j))*wanc(m2,j)
    enddo
  enddo
enddo
! compute S^{-1/2}
call isqrtzhe(nwantot,s,ierr)
ierr_=ierr
! compute Wannier function expansion coefficients
wanc_ort=zzero
if (ierr.eq.0) then
  do m1=1,nwantot
    do m2=1,nwantot
      wanc_ort(m1,:)=wanc_ort(m1,:)+wanc(m2,:)*s(m2,m1)
    enddo
  enddo
  wanc=wanc_ort
endif
deallocate(s,wanc_ort)
return
end subroutine

subroutine wan_gen_h(esv,wanc,h)
use modmain
implicit none
real(8), intent(in) :: esv(nstsv)
complex(8), intent(in) :: wanc(nwantot,nstsv)
complex(8), intent(out) :: h(nwantot,nwantot)
!
integer j,n1,n2
!
h(:,:)=zzero
do n1=1,nwantot
  do n2=1,nwantot
    do j=1,nstsv
      h(n1,n2)=h(n1,n2)+dconjg(wanc(n1,j))*wanc(n2,j)*esv(j)
    enddo
  enddo
enddo
return
end subroutine

! use two Fermi functions to make a smooth energy interval [e1,e2] with
!  the finite width w1,w2
real(8) function smoothstep(e,e1,e2,w1,w2)
implicit none
real(8), intent(in) :: e
real(8), intent(in) :: e1
real(8), intent(in) :: e2
real(8), intent(in) :: w1
real(8), intent(in) :: w2
smoothstep=1.d0/(exp((e-e2)/w2)+1.d0)+1.d0/(exp((-e+e1)/w1)+1.d0)-1.d0
if (smoothstep.lt.1d-16) smoothstep=0.d0 
return
end function


subroutine get_wanval(x,wanval,rmax,rcutoff)
use modmain
implicit none
real(8), intent(inout) :: x(3)
complex(8), intent(out) :: wanval(nspinor)
real(8), optional, intent(in) :: rmax
real(8), optional, intent(in) :: rcutoff
!
integer is,ia,ias,ir0,io,l,lm,ig,ispn
integer ntr(3),ik,ic
real(8) x0(3),vtc(3),vr0(3),r0,tp(2)
real(8) ur(0:lmaxapw,nufrmax),dr
complex(8) zt1,zt2,ylm(lmmaxapw)
logical, external :: vrinmt2
complex(8), external :: ylm_val
complex(8) expigr(m_ngvec)
complex(8) zm1(lmmaxapw,nufrmax,nspinor),zm2(nspinor)
!
wanval=zzero
x0(:)=x(:)-m_wanpos(:)
if (present(rmax)) then
  if (sum(x0(:)**2).gt.(rmax**2+1d-8)) then
    write(*,'("Error(s_get_wanval2): outside of rmax")')
    write(*,'("  rmax : ",G18.10)')rmax
    write(*,'("   x   : ",3G18.10)')x0
    write(*,'("  |x|  : ",G18.10)')sqrt(sum(x0(:)**2))
    call pstop
  endif
endif
if (present(rcutoff)) then
  if (sum(x0(:)**2).gt.(rcutoff**2)) return
endif

if (vrinmt2(x,is,ia,ntr,ir0,vr0,dr)) then
  ias=idxas(ia,is)
  ic=ias2ic(ias)
  call sphcrd(vr0,r0,tp)
  call genylm(lmaxapw,tp,ylm)
  vtc(:)=ntr(1)*avec(:,1)+ntr(2)*avec(:,2)+ntr(3)*avec(:,3)
  ur=0.d0
  do l=0,lmaxapw
    do io=1,nufr(l,is)
      ur(l,io)=ufr(ir0,l,io,ic)+dr*(ufr(ir0+1,l,io,ic)-ufr(ir0,l,io,ic))
    enddo !io
  enddo !l
  zm1=zzero
  do ik=1,nkptnr
    zt1=exp(zi*dot_product(vkcnr(:,ik),vtc(:)))*wkptnr(ik)
    do ispn=1,nspinor
      zm1(:,:,ispn)=zm1(:,:,ispn)+zt1*m_wann_unkmt(:,:,ias,ispn,ik)
    enddo
  enddo
  do ispn=1,nspinor
    do lm=1,lmmaxapw
      l=lm2l(lm)
      do io=1,nufr(l,is)
        wanval(ispn)=wanval(ispn)+zm1(lm,io,ispn)*ur(l,io)*ylm(lm)
      enddo !io
    enddo !lm
  enddo !ispn
else
  do ig=1,m_ngvec
    expigr(ig)=exp(zi*dot_product(x(:),vgc(:,ig)))  
  enddo
  do ik=1,nkptnr
    zt1=wkptnr(ik)*exp(zi*dot_product(vkcnr(:,ik),x(:)))/sqrt(omega)
    zm2=zzero
! TODO: check if this can be optimized by switching order of indexes
    do ispn=1,nspinor
      do ig=1,m_ngknr(ik)
      zt2=expigr(m_igkignr(ig,ik))
      !do ispn=1,nspinor
        zm2(ispn)=zm2(ispn)+zt2*m_wann_unkit(ig,ispn,ik)
      enddo
    enddo
    do ispn=1,nspinor
      wanval(ispn)=wanval(ispn)+zt1*zm2(ispn)
    enddo
  enddo
endif
return
end subroutine

end module
