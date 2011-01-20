subroutine test_spherical
use modmain
use mod_nrkp
use mod_wannier
use mod_madness
implicit none
integer i,ik,ikloc,ig
real(8) d1,d2,vrc(3)

integer ntp,it,itp,ias,ias1,ia,is,ir,lm,ir2,l
real(8), allocatable :: xmt_(:,:,:)
real(8), allocatable :: tp(:,:)
real(8), allocatable :: rlmb(:,:)
real(8), allocatable :: wtp(:)
real(8), allocatable :: spx(:,:)
real(8) t1,r1,r2,rr1,rr2
real(8), allocatable :: s_wantp(:,:)
real(8), allocatable :: s_wanlm(:,:)
real(8), allocatable :: s_pothlm(:,:)
integer s_nr
real(8) s_a,s_b,x1,x2
real(8), allocatable :: s_r(:),rlm(:)
integer lmaxwan
integer lmmaxwan

call init0
call init1
write(*,*)"Hello from Elk from process ",iproc
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
call getufr
call genufrp
wproc=.false.
call genwfnr(-1,.false.)

if (allocated(m_ngknr)) deallocate(m_ngknr)
allocate(m_ngknr(nkptnr))
m_ngknr=0
if (allocated(m_igkignr)) deallocate(m_igkignr)
allocate(m_igkignr(ngkmax,nkptnr))
m_igkignr=0
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  m_ngknr(ik)=ngknr(ikloc)
  m_igkignr(:,ik)=igkignr(:,ikloc)
enddo
call mpi_grid_reduce(m_ngknr(1),nkptnr,dims=(/dim_k/),all=.true.)
call mpi_grid_reduce(m_igkignr(1,1),ngkmax*nkptnr,dims=(/dim_k/),all=.true.)
m_ngvec=0
do ik=1,nkptnr
  do ig=1,m_ngknr(ik)
    m_ngvec=max(m_ngvec,m_igkignr(ig,ik))
  enddo
enddo

if (allocated(m_wann_unkmt)) deallocate(m_wann_unkmt)
allocate(m_wann_unkmt(lmmaxvr,nufrmax,natmtot,nspinor,nkptnr))
if (allocated(m_wann_unkit)) deallocate(m_wann_unkit)
allocate(m_wann_unkit(ngkmax,nspinor,nkptnr))







lmaxwan=10
lmmaxwan=(lmaxwan+1)**2

ntp=266
allocate(tp(2,ntp))
allocate(rlmb(ntp,lmmaxwan))
allocate(wtp(ntp))
allocate(rlm(lmmaxwan))
! Lebedev mesh
allocate(spx(3,ntp))
call leblaik(ntp,spx,wtp)                                                                   
do itp=1,ntp                   
  wtp(itp)=wtp(itp)*fourpi
  call sphcrd(spx(:,itp),t1,tp(:, itp))                                                     
enddo 
! generate spherical harmonics
do itp=1,ntp 
  call genrlm(lmaxwan,tp(1,itp),rlm)  
  do lm=1,lmmaxwan
    rlmb(itp,lm)=rlm(lm)*wtp(itp)
  enddo
enddo
!deallocate(wtp,tp)

s_nr=2000
x1=0.001d0
x2=12.d0
s_b=log(x2/x1)/(s_nr-1)
s_a=x1/exp(s_b)


allocate(s_wantp(ntp,s_nr))
allocate(s_wanlm(lmmaxwan,s_nr))
allocate(s_pothlm(lmmaxwan,s_nr))

allocate(s_r(s_nr))
do ir=1,s_nr
  s_r(ir)=s_a*exp(ir*s_b)
enddo


call elk_load_wann_unk(4)
s_wantp=0.d0
do ir=1,s_nr
  do itp=1,ntp
    vrc(:)=spx(:,itp)*s_r(ir)
    call elk_wan_rho(10.d0,vrc,s_wantp(itp,ir))
  enddo
enddo
call dgemm('T','N',lmmaxwan,s_nr,ntp,1.d0,rlmb,ntp,s_wantp,ntp,0.d0,&
  s_wanlm,lmmaxwan)

d1=0.d0
do ir=1,s_nr-1
  d1=d1+0.5d0*(s_wanlm(1,ir)*s_r(ir)**2+s_wanlm(1,ir+1)*s_r(ir+1)**2)*(s_r(ir+1)-s_r(ir))
enddo
d1=d1*fourpi*y00

write(*,*)"norm : ",d1

vrc=0.d0
do ir=1,s_nr
  vrc(1)=s_r(ir)
  call sphcrd(vrc,t1,tp(:,1))
  call genrlm(lmaxwan,tp(1,1),rlm)
  d1=0.d0
  do lm=1,lmmaxwan
    d1=d1+rlm(lm)*s_wanlm(lm,ir)
  enddo
  write(100,*)s_r(ir),d1
enddo


s_pothlm=0.d0
do lm=1,lmmaxwan
  l=lm2l(lm)
  do ir=1,s_nr
    r1=s_r(ir)
    t1=0.d0
    do ir2=1,s_nr-1
      r2=s_r(ir2)
      if (r1.lt.r2) then
        rr1=r1**l/r2**(l+1)
      else
        rr1=r2**l/r1**(l+1)
      endif
      r2=s_r(ir2+1)
      if (r1.lt.r2) then
        rr2=r1**l/r2**(l+1)
      else
        rr2=r2**l/r1**(l+1)
      endif
      t1=t1+0.5d0*(rr1*s_wanlm(lm,ir2)*s_r(ir2)**2+&
        rr2*s_wanlm(lm,ir2+1)*s_r(ir2+1)**2)*(s_r(ir2+1)-s_r(ir2))
    enddo !ir2
    s_pothlm(lm,ir)=t1*fourpi/(2*l+1)
  enddo !ir
enddo !lm
      
vrc=0.d0
do ir=1,s_nr
  vrc(1)=s_r(ir)
  call sphcrd(vrc,t1,tp(:,1))
  call genrlm(lmaxwan,tp(1,1),rlm)
  d1=0.d0
  do lm=1,lmmaxwan
    d1=d1+rlm(lm)*s_pothlm(lm,ir)
  enddo
  write(101,*)s_r(ir),d1
enddo

d1=0.d0
do lm=1,lmmaxwan
  t1=0.d0
  do ir=1,s_nr-1
    t1=t1+0.5d0*(s_wanlm(lm,ir)*s_pothlm(lm,ir)*s_r(ir)**2+&
      s_wanlm(lm,ir+1)*s_pothlm(lm,ir+1)*s_r(ir+1)**2)*(s_r(ir+1)-s_r(ir))
  enddo
  write(*,*)t1
  d1=d1+t1
enddo
write(*,*)"Hartree potential : ",d1


call bstop
return
end
