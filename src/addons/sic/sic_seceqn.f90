subroutine sic_seceqn(ikloc,evecfv,evecsv)
use modmain
use mod_lf
use mod_hdf5
! arguments
implicit none
integer, intent(in) :: ikloc
complex(8), intent(in) :: evecfv(nmatmax,nstfv,nspnfv)
complex(8), intent(inout) :: evecsv(nstsv,nstsv)
! local vars
logical exist
complex(8), allocatable :: hunif(:,:),z2(:,:)
complex(8), allocatable :: work(:)
real(8), allocatable :: rwork(:)
real(8), allocatable :: vn(:),hn(:)
integer i,n,ik,j,lwork,i1,i2,info
complex(8), allocatable :: wfsvmt(:,:,:,:,:)
complex(8), allocatable :: wfsvit(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfmt(:,:,:,:)
complex(8), allocatable :: wfir(:,:)
complex(8), allocatable :: wvp(:,:),u(:,:)
integer lm,ias,is,ispn,ig,ir,io
complex(8), allocatable :: hwank(:,:)
complex(8), allocatable :: vwank_1(:,:)
complex(8), allocatable :: vwank_2(:,:)
complex(8) expikt
integer v1l(3),n1,j1,i3
real(8) vtrc(3),v2(3),v3(3)
complex(8) zt1,expikr
real(8), parameter :: epsherm=1d-8

inquire(file="sic.hdf5",exist=exist)
if (.not.exist) return

allocate(hunif(nstsv,nstsv))
allocate(z2(nstsv,nstsv))
lwork=2*nstsv
allocate(rwork(3*nstsv))
allocate(work(lwork))

allocate(vn(nwann))
allocate(hn(nwann))
do i=1,nmegqwan
  if ((imegqwan(1,i).eq.imegqwan(2,i)).and.imegqwan(3,i).eq.0.and.&
    imegqwan(4,i).eq.0.and.imegqwan(5,i).eq.0) then
    vn(imegqwan(1,i))=dreal(vwan(i))
    !hn(imegqwan(1,i))=dreal(hwan(i))
  endif
enddo

ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)

allocate(wfsvmt(lmmaxvr,nufrmax,natmtot,nspinor,nstsv))
allocate(wfsvit(ngkmax,nspinor,nstsv))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wfmt(lmmaxvr,nrmtmax,natmtot,nspinor))
allocate(wfir(ngrtot,nspinor))
wfsvmt=zzero
wfsvit=zzero
! get apw coeffs 
call match(ngk(1,ik),gkc(:,1,ikloc),tpgkc(:,:,1,ikloc), &
 sfacgk(:,:,1,ikloc),apwalm)
! generate wave functions in muffin-tins
call genwfsvmt(lmaxvr,lmmaxvr,ngk(1,ik),evecfv,evecsv,apwalm,wfsvmt)
! generate wave functions in interstitial
call genwfsvit(ngk(1,ik),evecfv,evecsv,wfsvit)

! compute <W_n|V_n|\Psi_{jk}>
allocate(wvp(nwann,nstsv),u(nwann,nstsv))
wvp=zzero
u=zzero
do j=1,nstsv
  wfmt=zzero
  wfir=zzero
  do ispn=1,nspinor
    do ias=1,natmtot
      is=ias2is(ias)
      do ir=1,nrmt(is)
        do lm=1,lmmaxvr
          do io=1,nufr(lm2l(lm),is)
            wfmt(lm,ir,ias,ispn)=wfmt(lm,ir,ias,ispn)+&
              ufr(ir,lm2l(lm),io,ias2ic(ias))*wfsvmt(lm,io,ias,ispn,j)
          enddo
        enddo
      enddo
    enddo !ias
    do ig=1,ngk(1,ik)
      wfir(igfft(igkig(ig,1,ikloc)),ispn)=wfsvit(ig,ispn,j)
    enddo
    call zfftifc(3,ngrid,1,wfir(:,ispn))
    wfir(:,ispn)=wfir(:,ispn)/sqrt(omega)
  enddo !ispn
  ir=0
  do i3=0,ngrid(3)-1
    v2(3)=dble(i3)/dble(ngrid(3))
    do i2=0,ngrid(2)-1
      v2(2)=dble(i2)/dble(ngrid(2))
      do i1=0,ngrid(1)-1
        v2(1)=dble(i1)/dble(ngrid(1))
        ir=ir+1
        call r3mv(avec,v2,v3)
        expikr=exp(zi*dot_product(vkc(:,ik),v3(:)))
        wfir(ir,:)=expikr*wfir(ir,:)
      enddo
    enddo
  enddo
  do n=1,nwann
    do ispn=1,nspinor
      wvp(n,j)=wvp(n,j)+lf_dotblh(.true.,vkc(1,ik),vwanmt(1,1,1,1,ispn,n),&
        vwanir(1,1,ispn,n),wfmt(1,1,1,ispn),wfir(1,ispn))
    enddo
  enddo
enddo !j
deallocate(wfsvmt,wfsvit,apwalm,wfmt,wfir)
u=wann_c(:,:,ikloc)
 
! compute H_{nn'}(k) and V_{nn'}(k)
allocate(hwank(nwann,nwann))
allocate(vwank_1(nwann,nwann))
allocate(vwank_2(nwann,nwann))
hwank=zzero
vwank_1=zzero
vwank_2=zzero
do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  v1l(:)=imegqwan(3:5,i)
  vtrc(:)=v1l(1)*avec(:,1)+v1l(2)*avec(:,2)+v1l(3)*avec(:,3)
  expikt=exp(zi*dot_product(vkc(:,ik),vtrc(:)))
  !hwank(n,n1)=hwank(n,n1)+expikt*hwan(i)
  vwank_1(n,n1)=vwank_1(n,n1)+expikt*vwan(i)
  vwank_2(n1,n)=vwank_2(n1,n)+dconjg(expikt*vwan(i))
enddo

do n=1,nwann
  do n1=1,nwann
    do j=1,nstsv
      hwank(n,n1)=hwank(n,n1)+&
        dconjg(wann_c(n,j,ikloc))*wann_c(n1,j,ikloc)*evalsv(j,ik)
        !dconjg(wann_c(n1,j,ikloc))*wann_c(n,j,ikloc)*evalsv(j,ik)
    enddo
  enddo
enddo

! setup unified Hamiltonian
hunif=zzero
do j=1,nstsv
! 1-st term : diagonal H^{LDA} 
  hunif(j,j)=evalsv(j,ik) 
  do j1=1,nstsv
! 2-nd term : -\sum'_{\alpha,\alpha'} P_{\alpha} H^{LDA} P_{\alpha'}
    do n=1,nwann
      do n1=1,nwann
!        if (n.ne.n1) then
          hunif(j,j1)=hunif(j,j1)-hwank(n,n1)*u(n,j)*dconjg(u(n1,j1))
!        endif
      enddo
    enddo
! 3-rd term : \sum_{alpha} P_{\alpha} V_{\alpha} P_{\alpha})
    do n=1,nwann
      hunif(j,j1)=hunif(j,j1)+u(n,j)*dconjg(u(n,j1))*(vn(n)+wann_ene(n))
    enddo
! 4-th and 5-th terms : \sum_{\alpha} P_{\alpha} V_{\alpha} Q + 
!                       \sum_{\alpha} Q V_{\alpha} P_{\alpha} 
!   where Q=1-\sum_{\alpha'}P_{\alpha'}
    do n=1,nwann
      hunif(j,j1)=hunif(j,j1)+u(n,j)*wvp(n,j1)+dconjg(wvp(n,j)*u(n,j1))
      do n1=1,nwann
        hunif(j,j1)=hunif(j,j1)-vwank_1(n,n1)*u(n,j)*dconjg(u(n1,j1))-&
          vwank_2(n1,n)*u(n1,j)*dconjg(u(n,j1))
      enddo
    enddo
  enddo !j1
enddo !j

! check hermiticity
do j=1,nstsv
  do j1=1,nstsv
    if (abs(hunif(j,j1)-dconjg(hunif(j1,j))).gt.epsherm) then
      write(*,*)
      write(*,'("Warning(sic_seceqn) : unified Hamiltonian is not hermitian")')
    endif
  enddo
enddo

if (ndmag.eq.1) then
! collinear: block diagonalise H
  call zheev('V','U',nstfv,hunif(1,1),nstsv,evalsv(1,ik),work,lwork,rwork,info)
  if (info.ne.0) goto 20
  i=nstfv+1
  call zheev('V','U',nstfv,hunif(i,i),nstsv,evalsv(i,ik),work,lwork,rwork,info)
  if (info.ne.0) goto 20
  do i=1,nstfv
    do j=1,nstfv
      hunif(i,j+nstfv)=0.d0
      hunif(i+nstfv,j)=0.d0
    end do
  end do
else
! non-collinear or spin-unpolarised: full diagonalisation
  call zheev('V','U',nstsv,hunif,nstsv,evalsv(1,ik),work,lwork,rwork,info)
  if (info.ne.0) goto 20
endif
z2=zzero
do i1=1,nstsv
  do i2=1,nstsv
    do j=1,nstsv
      z2(i2,i1)=z2(i2,i1)+dconjg(hunif(j,i1))*evecsv(i2,j)
    enddo
  enddo
enddo
evecsv=z2

deallocate(hunif)
deallocate(z2)
deallocate(rwork)
deallocate(work)
deallocate(vn,hn)
deallocate(hwank,vwank_1,vwank_2)
deallocate(wvp,u)

call genwann(ikloc,evecfv,evecsv)
return
20 continue
write(*,*)
write(*,'("Error(sic_seceqn): diagonalisation of the second-variational &
 &Hamiltonian failed")')
write(*,'(" for k-point ",I8)') ik
write(*,'(" ZHEEV returned INFO = ",I8)') info
write(*,*)
call pstop
end