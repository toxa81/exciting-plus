subroutine sic_hunif(ikloc,evecfv,hunif)
use modmain
use mod_lf
use mod_hdf5
! arguments
implicit none
integer, intent(in) :: ikloc
complex(8), intent(in) :: evecfv(nmatmax,nstfv)
complex(8), intent(inout) :: hunif(nstsv,nstsv)
! local variables
logical exist
real(8), allocatable :: vn(:)
integer i,n,ik,i1,i2
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: wfmt(:,:,:)
complex(8), allocatable :: wfir(:)
complex(8), allocatable :: vwank(:,:)
complex(8), allocatable :: vwank_sym(:,:)
complex(8), allocatable :: a(:,:,:)
complex(8), allocatable :: b(:,:,:)
complex(8), allocatable :: zm1(:,:)
integer ias,ispn,ig,ir
complex(8) expikt
integer vtrl(3),n1,j1,j2,i3,ispn1,ispn2,jst1,jst2,ist,n2
real(8) vtrc(3),v2(3),v3(3)
complex(8) expikr
real(8) t1
real(8), parameter :: epsherm=1d-8
!character*100 fname

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
allocate(a(nwann,nstfv,nspinor))
allocate(b(nwann,nstfv,nspinor))

a=zzero
b=zzero
! get apw coeffs 
call match(ngk(1,ik),gkc(:,1,ikloc),tpgkc(:,:,1,ikloc),sfacgk(:,:,1,ikloc),&
  apwalm)
! compute a=<W_n|\phi_{jk}> and b=<W_n|V_n|\phi_{jk}> where phi(r) are firt-
!  variational Bloch wave-functions
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
      a(n,ist,ispn)=lf_dot_blh(.true.,vkc(1,ik),wanmt(1,1,1,1,ispn,n),&
        wanir(1,1,ispn,n),wfmt,wfir)
      b(n,ist,ispn)=lf_dot_blh(.true.,vkc(1,ik),wvmt(1,1,1,1,ispn,n),&
        wvir(1,1,ispn,n),wfmt,wfir)
    enddo
  enddo
enddo !ist
deallocate(apwalm,wfmt,wfir)
! compute V_{nn'}(k)
allocate(vwank(nwann,nwann))
allocate(vwank_sym(nwann,nwann))
vwank=zzero
do i=1,nmegqwan
  n1=imegqwan(1,i)
  n2=imegqwan(2,i)
  vtrl(:)=imegqwan(3:5,i)
  vtrc(:)=vtrl(1)*avec(:,1)+vtrl(2)*avec(:,2)+vtrl(3)*avec(:,3)
  expikt=exp(zi*dot_product(vkc(:,ik),vtrc(:)))
  vwank(n1,n2)=vwank(n1,n2)+expikt*vwanme(i)
enddo
! symmetrize the matrix
do n1=1,nwann
  do n2=1,nwann
    vwank_sym(n1,n2)=0.5d0*(vwank(n1,n2)+dconjg(vwank(n2,n1)))
  enddo
enddo
! compute H_{nn'}^{0}(k); remember that on entrance hunif=H0
allocate(zm1(nwann,nstsv))
call zgemm('N','N',nwann,nstsv,nstsv,zone,a,nwann,hunif,nstsv,zzero,zm1,nwann)
call zgemm('N','C',nwann,nwann,nstsv,zone,zm1,nwann,a,nwann,zzero,&
  sic_wann_h0k(1,1,ikloc),nwann)
deallocate(zm1)

!if (mpi_grid_root((/dim2/))) then
!  write(fname,'("h0_n",I2.2,"_k",I4.4".txt")')nproc,ik
!  call wrmtrx(fname,nstsv,nstsv,hunif,nstsv)
!  write(fname,'("sic_wann_h0k_n",I2.2,"_k",I4.4".txt")')nproc,ik
!  call wrmtrx(fname,nwann,nwann,sic_wann_h0k,nwann)
!endif

! setup unified Hamiltonian
! 1-st term: LDA Hamiltonian itself
do j1=1,nstfv
  do ispn1=1,nspinor
    jst1=j1+(ispn1-1)*nstfv
    do j2=1,nstfv
      do ispn2=1,nspinor
        jst2=j2+(ispn2-1)*nstfv
! 2-nd term : -\sum_{\alpha,\alpha'} P_{\alpha} H^{LDA} P_{\alpha'}
        do n1=1,nwann
          do n2=1,nwann
            hunif(jst1,jst2)=hunif(jst1,jst2)-sic_wann_h0k(n1,n2,ikloc)*&
              dconjg(a(n1,j1,ispn1))*a(n2,j2,ispn2)
          enddo
        enddo
! 3-rd term : \sum_{alpha} P_{\alpha} H^{LDA} P_{\alpha}
! 4-th term : \sum_{alpha} P_{\alpha} V_{\alpha} P_{\alpha}
        do n=1,nwann
          hunif(jst1,jst2)=hunif(jst1,jst2)+dconjg(a(n,j1,ispn1))*&
            a(n,j2,ispn2)*(vn(n)+sic_wann_e0(n))
        enddo
! 5-th term : \sum_{\alpha} P_{\alpha} V_{\alpha} Q + 
!             \sum_{\alpha} Q V_{\alpha} P_{\alpha} 
!  where Q=1-\sum_{\alpha'}P_{\alpha'}
        do n=1,nwann
          hunif(jst1,jst2)=hunif(jst1,jst2)+&
            dconjg(a(n,j1,ispn1))*b(n,j2,ispn2)+&
            dconjg(b(n,j1,ispn1))*a(n,j2,ispn2)
        enddo
        do n1=1,nwann
          do n2=1,nwann
            hunif(jst1,jst2)=hunif(jst1,jst2)-&
              vwank_sym(n1,n2)*dconjg(a(n1,j1,ispn1))*a(n2,j2,ispn2)-&
              dconjg(vwank_sym(n1,n2))*dconjg(a(n2,j1,ispn1))*a(n1,j2,ispn2)
          enddo
        enddo
      enddo !ispn2
    enddo !j2
  enddo !ispn1
enddo !j1

!if (mpi_grid_root((/dim2/))) then
!  write(fname,'("hunif_n",I2.2,"_k",I4.4".txt")')nproc,ik
!  call wrmtrx(fname,nstsv,nstsv,hunif,nstsv)
!endif

! check hermiticity
if (mpi_grid_root((/dim2/))) then  
  do j1=1,nstsv
    do j2=1,nstsv
      t1=abs(hunif(j1,j2)-dconjg(hunif(j2,j1)))
      if (t1.gt.epsherm) then
        write(*,*)
        write(*,'("Warning(sic_hunif) : unified Hamiltonian is not hermitian")')
        write(*,'("  k-point : ",I4)')ik
        write(*,'("  j1, j2, diff : ",2I4,G18.10)')j1,j2,t1
      endif
    enddo
  enddo
endif
deallocate(vn)
deallocate(a,b)
deallocate(vwank,vwank_sym)
return
end