complex(8) function inner_product(ntr,ntrloc,vtrl,t,f1mt,f1ir,f2mt,f2ir)
use modmain
implicit none
integer, intent(in) :: ntr
integer, intent(in) :: ntrloc
integer, intent(in) :: vtrl(3,ntr)
integer, intent(in) :: t(3)
complex(8), intent(in) :: f1mt(lmmaxvr,nrmtmax,natmtot,nspinor,ntrloc)
complex(8), intent(in) :: f1ir(ngrtot,nspinor,ntrloc)
complex(8), intent(in) :: f2mt(lmmaxvr,nrmtmax,natmtot,nspinor,ntrloc)
complex(8), intent(in) :: f2ir(ngrtot,nspinor,ntrloc)
complex(8) zprod
integer ispn,itrloc
complex(8), external :: zfinp_
! generates <f1_0|f2_t>
zprod=zzero
do itrloc=1,ntrloc
  do ispn=1,nspinor
    zprod=zprod+zfinp_(.true.,f1mt(1,1,1,ispn,itrloc),&
      f2mt(1,1,1,ispn,itrloc),f1ir(1,ispn,itrloc),f2ir(1,ispn,itrloc))
  enddo
enddo
call mpi_grid_reduce(zprod,dims=(/dim2/))
inner_product=zprod
return
end