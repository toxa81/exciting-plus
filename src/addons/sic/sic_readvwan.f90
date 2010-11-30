subroutine sic_readvwan
use modmain
use mod_sic
use mod_hdf5
use mod_linresp
use mod_wannier
implicit none
integer it,n,ispn
character*20 c1,c2,c3
character*100 path
logical exist
complex(8), allocatable :: fmt(:,:,:)
complex(8), allocatable :: fir(:)

! retutn if wannier functions and potential are initialized
if (tsic_wv) return
inquire(file="sic.hdf5",exist=exist)
if (.not.exist) return

call hdf5_read("sic.hdf5","/","nmegqwan",nmegqwan)
if (.not.allocated(imegqwan)) allocate(imegqwan(5,nmegqwan))
call hdf5_read("sic.hdf5","/","imegqwan",imegqwan(1,1),(/5,nmegqwan/))
if (allocated(vwanme)) deallocate(vwanme)
allocate(vwanme(nmegqwan))
call hdf5_read("sic.hdf5","/","vwanme",vwanme(1),(/nmegqwan/))
call hdf5_read("sic.hdf5","/","sic_etot_correction",sic_etot_correction)

allocate(fmt(lmmaxvr,nrmtmax,natmtot))
allocate(fir(ngrtot))
do n=1,nwantot
  do ispn=1,nspinor
    do it=1,ntr
      write(c1,'("n",I4.4)')n
      write(c2,'("t",I4.4)')it
      write(c3,'("s",I4.4)')ispn 
      path="/wann/"//trim(adjustl(c1))//"/"//trim(adjustl(c3))//"/"//&
        trim(adjustl(c2))
      if (mpi_grid_root()) then
        call hdf5_read("sic.hdf5",path,"wvmt",fmt(1,1,1),&
          (/lmmaxvr,nrmtmax,natmtot/))
        call hdf5_read("sic.hdf5",path,"wvir",fir(1),(/ngrtot/))
      endif
      call mpi_grid_bcast(fmt(1,1,1),lmmaxvr*nrmtmax*natmtot)
      call mpi_grid_bcast(fir(1),ngrtot)
      call sic_copy_mt_z(.true.,lmmaxvr,fmt,wvmt(1,1,it,ispn,n))
      call sic_copy_ir_z(.true.,fir,wvir(1,it,ispn,n))
      if (mpi_grid_root()) then
        call hdf5_read("sic.hdf5",path,"wanmt",fmt(1,1,1),&
          (/lmmaxvr,nrmtmax,natmtot/))
        call hdf5_read("sic.hdf5",path,"wanir",fir(1),(/ngrtot/))
      endif
      call mpi_grid_bcast(fmt(1,1,1),lmmaxvr*nrmtmax*natmtot)
      call mpi_grid_bcast(fir(1),ngrtot)
      call sic_copy_mt_z(.true.,lmmaxvr,fmt,wanmt(1,1,it,ispn,n))
      call sic_copy_ir_z(.true.,fir,wanir(1,it,ispn,n))
    enddo
  enddo
enddo
deallocate(fmt,fir)
tsic_wv=.true.
return
end
