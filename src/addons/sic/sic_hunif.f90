subroutine sic_hunif(ikloc,wann_ufv,evecfv,evecsv)
use modmain
use mod_lf
use mod_hdf5
! arguments
implicit none
integer, intent(in) :: ikloc
complex(8), intent(in) :: wann_ufv(nwann,nstfv,nspinor)
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(inout) :: evecsv(nstsv,nstsv)
! local variables
logical exist
complex(8), allocatable :: hunif(:,:)
complex(8), allocatable :: work(:)
real(8), allocatable :: rwork(:)
real(8), allocatable :: vn(:)
integer i,n,ik,j,lwork,i1,i2,info
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfmt(:,:,:)
complex(8), allocatable :: wfir(:)
complex(8), allocatable :: hwank(:,:)
complex(8), allocatable :: vwank(:,:)
complex(8), allocatable :: wanvblh(:,:,:)

integer ias,ispn,ig,ir
complex(8) expikt
integer vtrl(3),n1,j1,j2,i3,ispn1,ispn2,jst1,jst2,ist,n2
real(8) vtrc(3),v2(3),v3(3)
complex(8) expikr
real(8), parameter :: epsherm=1d-8

inquire(file="sic.hdf5",exist=exist)
if (.not.exist) return
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)

allocate(vn(nwann))
do i=1,nmegqwan
  if ((imegqwan(1,i).eq.imegqwan(2,i)).and.imegqwan(3,i).eq.0.and.&
    imegqwan(4,i).eq.0.and.imegqwan(5,i).eq.0) then
    vn(imegqwan(1,i))=dreal(vwanme(i))
  endif
enddo

allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(wfmt(lmmaxvr,nrmtmax,natmtot))
allocate(wfir(ngrtot))
allocate(wanvblh(nwann,nstfv,nspinor))
wanvblh=zzero
! get apw coeffs 
call match(ngk(1,ik),gkc(:,1,ikloc),tpgkc(:,:,1,ikloc),sfacgk(:,:,1,ikloc),&
  apwalm)
! compute <W_n|V_n|\phi_{jk}> where phi(r) is a firt-variational wave-function
do ist=1,nstfv
  wfmt=zzero
  wfir=zzero
  do ias=1,natmtot
    call wavefmt(1,lmaxvr,ias2is(ias),ias2ia(ias),ngk(1,ik),apwalm,&
      evecfv(:,ist),lmmaxvr,wfmt(1,1,ias))
  enddo
  do ig=1,ngk(1,ik)
    wfir(igfft(igkig(ig,1,ikloc)))=evecfv(ig,ist)
  enddo
  call zfftifc(3,ngrid,1,wfir)
  wfir(:)=wfir(:)/sqrt(omega)
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
        wfir(ir)=expikr*wfir(ir)
      enddo
    enddo
  enddo
  do n=1,nwann
    do ispn=1,nspinor
      wanvblh(n,ist,ispn)=lf_dot_lb(.true.,vkc(1,ik),wvmt(1,1,1,1,ispn,n),&
        wvir(1,1,ispn,n),wfmt,wfir)
    enddo
  enddo
enddo !ist
deallocate(apwalm,wfmt,wfir)
! compute H_{nn'}(k)
allocate(hwank(nwann,nwann))
hwank=zzero
do n1=1,nwann
  do n2=1,nwann
    do j=1,nstsv
      hwank(n1,n2)=hwank(n1,n2)+&
        dconjg(wann_c(n1,j,ikloc))*wann_c(n2,j,ikloc)*evalsv0(j,ik)
    enddo
  enddo
enddo
! compute V_{nn'}(k)
allocate(vwank(nwann,nwann))
vwank=zzero
do i=1,nmegqwan
  n1=imegqwan(1,i)
  n2=imegqwan(2,i)
  vtrl(:)=imegqwan(3:5,i)
  vtrc(:)=vtrl(1)*avec(:,1)+vtrl(2)*avec(:,2)+vtrl(3)*avec(:,3)
  expikt=exp(zi*dot_product(vkc(:,ik),vtrc(:)))
  vwank(n1,n2)=vwank(n1,n2)+expikt*vwanme(i)
enddo

!if (mpi_grid_root((/dim2/))) then
!  write(fname,'("wvp_n",I2.2,"_k",I4.4".txt")')nproc,ik
!  call wrmtrx(fname,nwann,nstsv,wvp,nwann)
!  write(fname,'("hwank_n",I2.2,"_k",I4.4".txt")')nproc,ik
!  call wrmtrx(fname,nwann,nwann,hwank,nwann)
!  write(fname,'("vwank_n",I2.2,"_k",I4.4".txt")')nproc,ik
!  call wrmtrx(fname,nwann,nwann,vwank,nwann)
!endif

! setup unified Hamiltonian
allocate(hunif(nstsv,nstsv))
! 1-st term: LDA Hamiltonian itself
hunif=hmltsv(:,:,ikloc)
do j1=1,nstfv
do ispn1=1,nspinor
  jst1=j1+(ispn1-1)*nstfv
  do j2=1,nstfv
  do ispn2=1,nspinor
    jst2=j2+(ispn2-1)*nstfv
! 2-nd term : -\sum'_{\alpha,\alpha'} P_{\alpha} H^{LDA} P_{\alpha'} = 
!   -\sum_{\alpha,\alpha'} P_{\alpha} H^{LDA} P_{\alpha'} +   <-- take this
!   +\sum_{\alpha} P_{\alpha} H^{LDA} P_{\alpha}
    do n1=1,nwann
      do n2=1,nwann
        hunif(jst1,jst2)=hunif(jst1,jst2)-hwank(n1,n2)*wann_ufv(n1,j1,ispn1)*&
          dconjg(wann_ufv(n2,j2,ispn2))
      enddo
    enddo
! 3-rd term : \sum_{alpha} P_{\alpha} V_{\alpha} P_{\alpha})
!  plus 4-th term: diagonal energy matrix element E_n
    do n=1,nwann
      hunif(jst1,jst2)=hunif(jst1,jst2)+wann_ufv(n,j1,ispn1)*&
        dconjg(wann_ufv(n,j2,ispn2))*(vn(n)+sic_wann_ene(n))
    enddo
! 5-th term : \sum_{\alpha} P_{\alpha} V_{\alpha} Q + 
!                       \sum_{\alpha} Q V_{\alpha} P_{\alpha} 
!   where Q=1-\sum_{\alpha'}P_{\alpha'}
    do n=1,nwann
      hunif(jst1,jst2)=hunif(jst1,jst2)+&
        wann_ufv(n,j1,ispn1)*wanvblh(n,j2,ispn2)+&
        dconjg(wanvblh(n,j1,ispn1)*wann_ufv(n,j2,ispn2))
    enddo
    do n1=1,nwann
      do n2=1,nwann
        hunif(jst1,jst2)=hunif(jst1,jst2)-&
          vwank(n1,n2)*wann_ufv(n1,j1,ispn1)*dconjg(wann_ufv(n2,j2,ispn2))-&
          dconjg(vwank(n1,n2))*wann_ufv(n2,j1,ispn1)*dconjg(wann_ufv(n1,j2,ispn2))
      enddo
    enddo
  enddo !ispn2
  enddo !j2
enddo !ispn1
enddo !j1

if (mpi_grid_root((/dim2/))) then! check hermiticity
  do j1=1,nstsv
    do j2=1,nstsv
      if (abs(hunif(j1,j2)-dconjg(hunif(j2,j1))).gt.epsherm) then
        write(*,*)
        write(*,'("Warning(sic_seceqn) : unified Hamiltonian is not hermitian")')
      endif
    enddo
  enddo
endif

if (mpi_grid_root((/dim2/))) then
  lwork=2*nstsv
  allocate(rwork(3*nstsv))
  allocate(work(lwork))
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
  evecsv=hunif
  deallocate(rwork)
  deallocate(work)
endif
call mpi_grid_bcast(evecsv(1,1),nstsv*nstsv,dims=(/dim2/))
call mpi_grid_bcast(evalsv(1,ik),nstsv,dims=(/dim2/))
deallocate(hunif)
deallocate(vn)
deallocate(wanvblh)
deallocate(hwank)
deallocate(vwank)
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