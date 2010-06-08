subroutine seceqn_sic(ikloc,evecfv,evecsv)
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
complex(8), allocatable :: z1(:,:),z2(:,:),z5(:,:)
complex(8), allocatable :: work(:)
real(8), allocatable :: rwork(:)
complex(8), allocatable :: vsic(:)
complex(8), allocatable :: h0wan(:)
real(8), allocatable :: vn(:),hnn(:)
integer i,n,ik,j,lwork,i1,i2,info
complex(8), allocatable :: wfsvmt(:,:,:,:,:)
complex(8), allocatable :: wfsvit(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfmt(:,:,:,:)
complex(8), allocatable :: wfir(:,:)
complex(8), allocatable :: a1(:,:),a2(:,:)
integer lm,ias,is,ispn,ig,ir,io
complex(8), allocatable :: hwan_k(:,:)
complex(8), allocatable :: vwan_k(:,:)
complex(8), allocatable :: vwan_2_k(:,:)
complex(8) expikt
integer v1l(3),n1,j1,i3
real(8) vtrc(3),v2(3),v3(3)
complex(8) zt1,expikr
complex(8), external :: zfinp_
real(8), parameter :: epsherm=1d-8


inquire(file="sic.hdf5",exist=exist)
if (.not.exist) return

allocate(z1(nstsv,nstsv))
allocate(z2(nstsv,nstsv))
lwork=2*nstsv
allocate(rwork(3*nstsv))
allocate(work(lwork))

call hdf5_read("sic.hdf5","/","nmegqwan",nmegqwan)
if (.not.allocated(imegqwan)) allocate(imegqwan(5,nmegqwan))
call hdf5_read("sic.hdf5","/","imegqwan",imegqwan(1,1),(/5,nmegqwan/))
allocate(vsic(nmegqwan))
allocate(h0wan(nmegqwan))
call hdf5_read("sic.hdf5","/","vsic",vsic(1),(/nmegqwan/))
call hdf5_read("sic.hdf5","/","h0wan",h0wan(1),(/nmegqwan/))
allocate(vn(nwann))
allocate(hnn(nwann))
do i=1,nmegqwan
  if ((imegqwan(1,i).eq.imegqwan(2,i)).and.imegqwan(3,i).eq.0.and.&
    imegqwan(4,i).eq.0.and.imegqwan(5,i).eq.0) then
    vn(imegqwan(1,i))=dreal(vsic(i))
    hnn(imegqwan(1,i))=dreal(h0wan(i))
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

allocate(a1(nwann,nstsv),a2(nwann,nstsv))
a1=zzero
a2=zzero
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
              ufr(ir,lm2l(lm),io,ias)*wfsvmt(lm,io,ias,ispn,j)
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
      a1(n,j)=a1(n,j)+lf_dotblh(.true.,vkc(1,ik),vwanmt(1,1,1,1,ispn,n),&
        vwanir(1,1,ispn,n),wfmt(1,1,1,ispn),wfir(1,ispn))
!      a2(n,j)=a2(n,j)+lf_dotblh(.false.,vkc(1,ik),wanmt(1,1,1,1,ispn,n),&
!        wanir(1,1,ispn,n),wfmt(1,1,1,ispn),wfir(1,ispn))
    enddo
  enddo
enddo !j
deallocate(wfsvmt,wfsvit,apwalm,wfmt,wfir)
a2=wann_c(:,:,ikloc)
 
allocate(hwan_k(nwann,nwann))
allocate(vwan_k(nwann,nwann))
allocate(vwan_2_k(nwann,nwann))
hwan_k=zzero
vwan_k=zzero
vwan_2_k=zzero


do i=1,nmegqwan
  n=imegqwan(1,i)
  n1=imegqwan(2,i)
  v1l(:)=imegqwan(3:5,i)
  vtrc(:)=v1l(1)*avec(:,1)+v1l(2)*avec(:,2)+v1l(3)*avec(:,3)
  expikt=exp(zi*dot_product(vkc(:,ik),vtrc(:)))
  hwan_k(n,n1)=hwan_k(n,n1)+expikt*h0wan(i)
  vwan_k(n,n1)=vwan_k(n,n1)+expikt*vsic(i)
  vwan_2_k(n1,n)=vwan_2_k(n1,n)+dconjg(expikt*vsic(i))
enddo
!do n=1,nwann
!  do n1=1,nwann
!    if (abs(hwan_k(n,n1)-dconjg(hwan_k(n1,n))).gt.epsherm) then
!      write(*,*)
!      write(*,'("Warning : hwan_k is not hermitian")')
!    endif
!    if (abs(vwan_k(n,n1)-dconjg(vwan_k(n1,n))).gt.epsherm) then
!      write(*,*)
!      write(*,'("Warning : vwan_k is not hermitian")')
!      write(*,'(" difference : ",G18.10)')abs(vwan_k(n,n1)-dconjg(vwan_k(n1,n)))
!    endif
!  enddo
!enddo

z1=zzero
do i=1,nstsv
  z1(i,i)=evalsv(i,ik)
enddo

do j=1,nstsv
  do j1=1,nstsv
    do n=1,nwann
      do n1=1,nwann
        z1(j,j1)=z1(j,j1)-hwan_k(n,n1)*a2(n,j)*dconjg(a2(n1,j1))
      enddo
    enddo
  enddo
enddo

do n=1,nwann
  do j=1,nstsv
    do j1=1,nstsv
      z1(j,j1)=z1(j,j1)+a2(n,j)*dconjg(a2(n,j1))*(vn(n)+hnn(n))
    enddo
  enddo
enddo

allocate(z5(nstsv,nstsv))
z5=zzero
do j=1,nstsv
  do j1=1,nstsv
    do n=1,nwann
      z5(j,j1)=z5(j,j1)+a2(n,j)*a1(n,j1)+dconjg(a1(n,j)*a2(n,j1))
    enddo
  enddo
enddo
do j=1,nstsv
  do j1=1,nstsv
    do n=1,nwann
      do n1=1,nwann
        z5(j,j1)=z5(j,j1)-vwan_k(n,n1)*a2(n,j)*dconjg(a2(n1,j1))-&
        vwan_2_k(n1,n)*a2(n1,j)*dconjg(a2(n,j1))
      enddo
    enddo
  enddo
enddo  

do j=1,nstsv
  do j1=1,nstsv
    z1(j,j1)=z1(j,j1)+z5(j,j1)
  enddo
enddo

do j=1,nstsv
  do j1=1,nstsv
    if (abs(z1(j,j1)-dconjg(z1(j1,j))).gt.epsherm) then
      write(*,*)
      write(*,'("Warning : unified Hamiltonian is not hermitian")')
    endif
  enddo
enddo

if (ndmag.eq.1) then
! collinear: block diagonalise H
  call zheev('V','U',nstfv,z1(1,1),nstsv,evalsv(1,ik),work,lwork,rwork,info)
  if (info.ne.0) goto 20
  i=nstfv+1
  call zheev('V','U',nstfv,z1(i,i),nstsv,evalsv(i,ik),work,lwork,rwork,info)
  if (info.ne.0) goto 20
  do i=1,nstfv
    do j=1,nstfv
      z1(i,j+nstfv)=0.d0
      z1(i+nstfv,j)=0.d0
    end do
  end do
else
! non-collinear or spin-unpolarised: full diagonalisation
  call zheev('V','U',nstsv,z1,nstsv,evalsv(1,ik),work,lwork,rwork,info)
  if (info.ne.0) goto 20
endif
z2=zzero
do i1=1,nstsv
  do i2=1,nstsv
    do j=1,nstsv
      z2(i2,i1)=z2(i2,i1)+dconjg(z1(j,i1))*evecsv(i2,j)
    enddo
  enddo
enddo
evecsv=z2

deallocate(z1)
deallocate(z2,z5)
deallocate(rwork)
deallocate(work)
deallocate(vsic,h0wan,vn,hnn)
deallocate(hwan_k,vwan_k,vwan_2_k)
deallocate(a1,a2)

call genwann(ikloc,evecfv,evecsv)
return
20 continue
write(*,*)
write(*,'("Error(seceqn_sic): diagonalisation of the second-variational &
 &Hamiltonian failed")')
write(*,'(" for k-point ",I8)') ik
write(*,'(" ZHEEV returned INFO = ",I8)') info
write(*,*)
call pstop
end