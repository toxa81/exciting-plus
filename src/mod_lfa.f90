! 
! Lattice Function Algebra (lfa) 
!
module mod_lfa
use modmain

! maximum lattice translation
integer trmax
! total number of translations
integer ntr
! translation vectors in lattice coordinates
integer, allocatable :: vtl(:,:)
! vector -> index map
integer, allocatable :: ivtit(:,:,:)

contains

subroutine lfa_init(trmax_)
implicit none
integer, intent(in) :: trmax_ 
integer i1,i2,i3,n
trmax=trmax_
ntr=(2*trmax+1)**3
if (allocated(vtl)) deallocate(vtl)
allocate(vtl(3,ntr))
if (allocated(ivtit)) deallocate(ivtit)
allocate(ivtit(-trmax:trmax,-trmax:trmax,-trmax:trmax))
n=0
do i1=-trmax,trmax
  do i2=-trmax,trmax
    do i3=-trmax,trmax
      n=n+1
      vtl(:,n)=(/i1,i2,i3/)
      ivtit(i1,i2,i3)=n
    enddo
  enddo
enddo
return
end subroutine

complex(8) function lfa_dotp(tsh,t,f1mt,f1ir,f2mt,f2ir)
use modmain
implicit none
logical, intent(in) :: tsh
integer, intent(in) :: t(3)
complex(8), intent(in) :: f1mt(lmmaxvr,nrmtmax,natmtot,*)
complex(8), intent(in) :: f1ir(ngrtot,*)
complex(8), intent(in) :: f2mt(lmmaxvr,nrmtmax,natmtot,*)
complex(8), intent(in) :: f2ir(ngrtot,*)
complex(8), allocatable :: f2mt_tmp(:,:,:)
complex(8), allocatable :: f2ir_tmp(:)

complex(8) zprod
integer ispn,itloc,it,jt,v1(3),v2(3),j,ntstep,ntrloc,itstep,ntloc1
integer jtloc,i,tag
logical l1
complex(8), external :: zfinp_
! generates <f1_0|f2_t>

zprod=zzero
allocate(f2mt_tmp(lmmaxvr,nrmtmax,natmtot))
allocate(f2ir_tmp(ngrtot))

j=0
ntstep=mpi_grid_map(ntr,dim2,x=j)
ntrloc=mpi_grid_map(ntr,dim2)
do itstep=1,ntstep
  f2mt_tmp=zzero
  f2ir_tmp=zzero
  do i=0,mpi_grid_size(dim2)-1
    ntloc1=mpi_grid_map(ntr,dim2,x=i)
    if (itstep.le.ntloc1) then
      it=mpi_grid_map(ntr,dim2,x=i,loc=itstep)
      v1(:)=vtl(:,it)
      v2(:)=v1(:)-t(:)
      l1=.false.
      if (v2(1).ge.-trmax.and.v2(1).le.trmax.and.&
          v2(2).ge.-trmax.and.v2(2).le.trmax.and.&
          v2(3).ge.-trmax.and.v2(3).le.trmax) then
        jt=ivtit(v2(1),v2(2),v2(3))
        l1=.true.
        jtloc=mpi_grid_map(ntr,dim2,glob=jt,x=j)
      endif
      if (l1.and.mpi_grid_x(dim2).eq.j.and.mpi_grid_x(dim2).ne.i) then
        tag=(itstep*mpi_grid_size(dim2)+i)*10
        call mpi_grid_send(f2mt(1,1,1,jtloc),lmmaxvr*nrmtmax*natmtot,&
          (/dim2/),(/i/),tag)
        call mpi_grid_send(f2ir(1,jtloc),ngrtot,(/dim2/),(/i/),tag+1)
      endif
      if (l1.and.mpi_grid_x(dim2).eq.i) then
        if (j.ne.i) then
          tag=(itstep*mpi_grid_size(dim2)+i)*10
          call mpi_grid_recieve(f2mt_tmp(1,1,1),lmmaxvr*nrmtmax*natmtot,&
            (/dim2/),(/j/),tag)
          call mpi_grid_recieve(f2ir_tmp(1),ngrtot,(/dim2/),(/j/),tag+1)
        else
          f2mt_tmp(:,:,:)=f2mt(:,:,:,jtloc)
          f2ir_tmp(:)=f2ir(:,jtloc)
        endif
      endif
    endif
  enddo !
  if (itstep.le.ntrloc) then
    zprod=zprod+zfinp_(tsh,f1mt(1,1,1,itstep),f2mt_tmp,f1ir(1,itstep),&
      f2ir_tmp)
  endif
enddo
deallocate(f2mt_tmp,f2ir_tmp)
call mpi_grid_reduce(zprod,dims=(/dim2/))
lfa_dotp=zprod
return
end function







end module