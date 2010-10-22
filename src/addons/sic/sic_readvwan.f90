subroutine sic_readvwan
use modmain
use mod_lf
use mod_hdf5
implicit none
integer itloc,it,n,ispn,nloc
character*20 c1,c2,c3
character*100 path
logical exist

inquire(file="sic.hdf5",exist=exist)
if (.not.exist) return

call hdf5_read("sic.hdf5","/","nmegqwan",nmegqwan)
if (.not.allocated(imegqwan)) allocate(imegqwan(5,nmegqwan))
call hdf5_read("sic.hdf5","/","imegqwan",imegqwan(1,1),(/5,nmegqwan/))
if (allocated(vwanme)) deallocate(vwanme)
allocate(vwanme(nmegqwan))
call hdf5_read("sic.hdf5","/","vwanme",vwanme(1),(/nmegqwan/))
call hdf5_read("sic.hdf5","/","sic_etot_correction",sic_etot_correction)

do itloc=1,ntrloc
  it=mpi_grid_map(ntr,dim_t,loc=itloc)
  do ispn=1,nspinor
    do nloc=1,nwannloc
      n=mpi_grid_map(nwann,dim_k,loc=nloc)
      write(c1,'("n",I4.4)')n
      write(c2,'("t",I4.4)')it
      write(c3,'("s",I4.4)')ispn 
      path="/wann/"//trim(adjustl(c1))//"/"//trim(adjustl(c3))//"/"//&
        trim(adjustl(c2))        
     call hdf5_read("sic.hdf5",path,"wvmt",wvmt(1,1,1,itloc,ispn,nloc),&
        (/lmmaxvr,nrmtmax,natmtot/))
      call hdf5_read("sic.hdf5",path,"wvir",wvir(1,itloc,ispn,nloc),&
        (/ngrtot/))
     call hdf5_read("sic.hdf5",path,"wanmt",wanmt(1,1,1,itloc,ispn,nloc),&
        (/lmmaxvr,nrmtmax,natmtot/))
      call hdf5_read("sic.hdf5",path,"wanir",wanir(1,itloc,ispn,nloc),&
        (/ngrtot/))
    enddo
  enddo
enddo
return
end