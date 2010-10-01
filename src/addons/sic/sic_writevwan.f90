subroutine sic_writevwan
use modmain
use mod_lf
use mod_hdf5
implicit none
integer n,ispn,it,i,itloc
character*12 c1,c2,c3
character*100 path

if (wproc) then
  call hdf5_create_file("sic.hdf5")
  call hdf5_create_group("sic.hdf5","/","wann")
  do n=1,nwann
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
  call hdf5_write("sic.hdf5","/","vwan",vwanme(1),(/nmegqwan/))
  call hdf5_write("sic.hdf5","/","sic_etot_correction",sic_etot_correction)
  call hdf5_write("sic.hdf5","/","wann_ene",wann_ene(1),(/nwann/))
endif
if (mpi_grid_side(dims=(/dim_t/))) then
  do i=0,mpi_grid_size(dim_t)-1
    if (mpi_grid_x(dim_t).eq.i) then
      do itloc=1,ntrloc
        it=mpi_grid_map(ntr,dim_t,loc=itloc)
        do ispn=1,nspinor
          do n=1,nwann
            write(c1,'("n",I4.4)')n
            write(c2,'("s",I4.4)')ispn
            write(c3,'("t",I4.4)')it
            path="/wann/"//trim(adjustl(c1))//"/"//trim(adjustl(c2))//"/"//&
              trim(adjustl(c3))       
            call hdf5_write("sic.hdf5",path,"wvmt",&
              wvmt(1,1,1,itloc,ispn,n),(/lmmaxvr,nrmtmax,natmtot/))
            call hdf5_write("sic.hdf5",path,"wvir",&
              wvir(1,itloc,ispn,n),(/ngrtot/))
          enddo
        enddo
      enddo
    endif
    call mpi_grid_barrier(dims=(/dim_t/))
  enddo
endif
return
end