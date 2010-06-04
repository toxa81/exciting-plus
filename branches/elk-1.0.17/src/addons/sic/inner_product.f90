complex(8) function inner_product(ntr,ntrloc,trmax,vtrl,ivtit,t,f1mt,f1ir,&
  f2mt,f2ir)
use modmain
implicit none
integer, intent(in) :: ntr
integer, intent(in) :: ntrloc
integer, intent(in) :: trmax
integer, intent(in) :: vtrl(3,ntr)
integer, intent(in) :: ivtit(-trmax:trmax,-trmax:trmax,-trmax:trmax)
integer, intent(in) :: t(3)
complex(8), intent(in) :: f1mt(lmmaxvr,nrmtmax,natmtot,nspinor,ntrloc)
complex(8), intent(in) :: f1ir(ngrtot,nspinor,ntrloc)
complex(8), intent(in) :: f2mt(lmmaxvr,nrmtmax,natmtot,nspinor,ntrloc)
complex(8), intent(in) :: f2ir(ngrtot,nspinor,ntrloc)
complex(8), allocatable :: f2mt_tmp(:,:,:,:)
complex(8), allocatable :: f2ir_tmp(:,:)

complex(8) zprod
integer ispn,itloc,it,jt,v1(3),v2(3),j,ntstep,itstep,ntloc1
integer jtloc,i,tag
logical l1
complex(8), external :: zfinp_
! generates <f1_0|f2_t>
zprod=zzero
allocate(f2mt_tmp(lmmaxvr,nrmtmax,natmtot,nspinor))
allocate(f2ir_tmp(ngrtot,nspinor))

j=0
ntstep=mpi_grid_map(ntr,dim2,x=j)
do itstep=1,ntstep
  f2mt_tmp=zzero
  f2ir_tmp=zzero
  do i=0,mpi_grid_size(dim2)-1
    ntloc1=mpi_grid_map(ntr,dim2,x=i)
    if (itstep.le.ntloc1) then
      it=mpi_grid_map(ntr,dim2,x=i,loc=itstep)
      v1(:)=vtrl(:,it)
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
        call mpi_grid_send(f2mt(1,1,1,1,jtloc),&
          lmmaxvr*nrmtmax*natmtot*nspinor,(/dim2/),(/i/),tag)
        call mpi_grid_send(f2ir(1,1,jtloc),&
          ngrtot*nspinor,(/dim2/),(/i/),tag+1)
      endif
      if (l1.and.mpi_grid_x(dim2).eq.i) then
        if (j.ne.i) then
          tag=(itstep*mpi_grid_size(dim2)+i)*10
          call mpi_grid_recieve(f2mt_tmp(1,1,1,1),&
            lmmaxvr*nrmtmax*natmtot*nspinor,(/dim2/),(/j/),tag)
          call mpi_grid_recieve(f2ir_tmp(1,1),&
            ngrtot*nspinor,(/dim2/),(/j/),tag+1)
        else
          f2mt_tmp(:,:,:,:)=f2mt(:,:,:,:,jtloc)
          f2ir_tmp(:,:)=f2ir(:,:,jtloc)
        endif
      endif
    endif
  enddo !
  if (itstep.le.ntrloc) then
    do ispn=1,nspinor
      zprod=zprod+zfinp_(.true.,f1mt(1,1,1,ispn,itstep),&
        f2mt_tmp(1,1,1,ispn),f1ir(1,ispn,itstep),f2ir_tmp(1,ispn))
    enddo
  endif
enddo
      
!do itloc=1,ntrloc
!  it=mpi_grid_map(ntr,dim2,loc=itloc)
!  v1(:)=vtrl(:,it)
!  write(*,*)'it=',it,'vt=',v1
!  v2(:)=v1(:)-t(:)
!  if (v2(1).ge.-trmax.and.v2(1).le.trmax.and.&
!      v2(2).ge.-trmax.and.v2(2).le.trmax.and.&
!      v2(3).ge.-trmax.and.v2(3).le.trmax) then
!    jt=ivtit(v2(1),v2(2),v2(3))
!    write(*,*)'  jt=',jt,'vt=',v2
!  endif  
!  
!  
!  !do ispn=1,nspinor
!  !  zprod=zprod+zfinp_(.true.,f1mt(1,1,1,ispn,itloc),&
!  !    f2mt(1,1,1,ispn,itloc),f1ir(1,ispn,itloc),f2ir(1,ispn,itloc))
!  !enddo
!enddo
deallocate(f2mt_tmp,f2ir_tmp)
call mpi_grid_reduce(zprod,dims=(/dim2/))
inner_product=zprod
return
end