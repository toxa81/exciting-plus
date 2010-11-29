subroutine sic_writevwan
use modmain
use mod_sic
use mod_hdf5
implicit none
integer n,ispn,it,i,itloc,nloc,j
character*12 c1,c2,c3
character*100 path
complex(8), allocatable :: fmt(:,:,:)
complex(8), allocatable :: fir(:)

allocate(fmt(lmmaxvr,nrmtmax,natmtot))
allocate(fir(ngrtot))

if (wproc) then
  call hdf5_create_file("sic.hdf5")
  call hdf5_create_group("sic.hdf5","/","wann")
  do n=1,nwantot
    path="/wann"
    write(c1,'("n",I4.4)')n
    call hdf5_create_group("sic.hdf5",path,trim(adjustl(c1)))   
    do ispn=1,nspinor
      path="/wann/"//trim(adjustl(c1))
      write(c2,'("s",I4.4)')ispn
      call hdf5_create_group("sic.hdf5",path,trim(adjustl(c2)))   
      do it=1,ntr
        path="/wann/"//trim(adjustl(c1))//"/"//trim(adjustl(c2))      
        write(c3,'("t",I4.4)')it
        call hdf5_create_group("sic.hdf5",path,trim(adjustl(c3)))
      enddo
    enddo
  enddo
  call hdf5_write("sic.hdf5","/","nmegqwan",nmegqwan)
  call hdf5_write("sic.hdf5","/","imegqwan",imegqwan(1,1),(/5,nmegqwan/))
  call hdf5_write("sic.hdf5","/","vwanme",vwanme(1),(/nmegqwan/))
  call hdf5_write("sic.hdf5","/","sic_etot_correction",sic_etot_correction)
endif
do n=1,nwantot
  do ispn=1,nspinor
    do it=1,ntr
      fmt=zzero
      fir=zzero
      call sic_copy_mt_z(.false.,lmmaxvr,fmt,wanmt(1,1,it,ispn,n))
      call sic_copy_ir_z(.false.,fir,wanir(1,it,ispn,n))
      call mpi_grid_reduce(fmt(1,1,1),lmmaxvr*nrmtmax*natmtot)
      call mpi_grid_reduce(fir(1),ngrtot)
      if (mpi_grid_root()) then
        write(c1,'("n",I4.4)')n
        write(c2,'("s",I4.4)')ispn
        write(c3,'("t",I4.4)')it
        path="/wann/"//trim(adjustl(c1))//"/"//trim(adjustl(c2))//"/"//&
          trim(adjustl(c3))       
        call hdf5_write("sic.hdf5",path,"wanmt",fmt(1,1,1),&
          (/lmmaxvr,nrmtmax,natmtot/))
        call hdf5_write("sic.hdf5",path,"wanir",fir(1),(/ngrtot/))
      endif
      fmt=zzero
      fir=zzero
      call sic_copy_mt_z(.false.,lmmaxvr,fmt,wvmt(1,1,it,ispn,n))
      call sic_copy_ir_z(.false.,fir,wvir(1,it,ispn,n))
      call mpi_grid_reduce(fmt(1,1,1),lmmaxvr*nrmtmax*natmtot)
      call mpi_grid_reduce(fir(1),ngrtot)
      if (mpi_grid_root()) then
        call hdf5_write("sic.hdf5",path,"wvmt",fmt(1,1,1),&
          (/lmmaxvr,nrmtmax,natmtot/))
        call hdf5_write("sic.hdf5",path,"wvir",fir(1),(/ngrtot/))
      endif
    enddo
  enddo
enddo
deallocate(fmt,fir)
return
end
