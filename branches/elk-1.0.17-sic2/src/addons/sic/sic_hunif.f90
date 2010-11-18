subroutine sic_hunif(ikloc,hunif)
use modmain
use mod_lf
use mod_hdf5
! arguments
implicit none
integer, intent(in) :: ikloc
complex(8), intent(inout) :: hunif(nstsv,nstsv)
! local variables
logical exist
integer i,n,ik,vtrl(3),n1,j1,j2,ispn1,ispn2,jst1,jst2,n2
real(8) vtrc(3),t1
real(8), allocatable :: vn(:)
complex(8) expikt
complex(8), allocatable :: vwank(:,:)
!complex(8), allocatable :: vwank_sym(:,:)
complex(8), allocatable :: a(:,:,:)
complex(8), allocatable :: b(:,:,:)
complex(8), allocatable :: zm1(:,:)
real(8), parameter :: epsherm=1d-10
!character*100 fname

inquire(file="sic.hdf5",exist=exist)
if (.not.exist) return
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)

! restore full hermitian matrix
do j1=2,nstsv
  do j2=1,j1-1
    hunif(j1,j2)=dconjg(hunif(j2,j1))
  enddo
enddo

allocate(vn(nwann))
do i=1,nmegqwan
  if ((imegqwan(1,i).eq.imegqwan(2,i)).and.imegqwan(3,i).eq.0.and.&
    imegqwan(4,i).eq.0.and.imegqwan(5,i).eq.0) then
    vn(imegqwan(1,i))=dreal(vwanme(i))
  endif
enddo

allocate(a(nwann,nstfv,nspinor))
allocate(b(nwann,nstfv,nspinor))
a(:,:,:)=sic_wb(:,:,:,ikloc)
b(:,:,:)=sic_wvb(:,:,:,ikloc)

! compute V_{nn'}(k)
allocate(vwank(nwann,nwann))
!allocate(vwank_sym(nwann,nwann))
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
!do n1=1,nwann
!  do n2=1,nwann
!    vwank_sym(n1,n2)=0.5d0*(vwank(n1,n2)+dconjg(vwank(n2,n1)))
!  enddo
!enddo
! compute H_{nn'}^{0}(k); remember that on input hunif=H0
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
              vwank(n1,n2)*dconjg(a(n1,j1,ispn1))*a(n2,j2,ispn2)-&
              dconjg(vwank(n1,n2))*dconjg(a(n2,j1,ispn1))*a(n1,j2,ispn2)
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
deallocate(vwank)
return
end