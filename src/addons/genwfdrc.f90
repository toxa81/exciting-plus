subroutine genwfdrc
use modmain
use mod_nrkp
use mod_addons_q
use mod_hdf5
implicit none
character*100 fname
character*20 kname,gname
integer ikloc,ik,i

call init0
call init1
if (.not.mpi_grid_in()) return
call readstate
call linengy
call genapwfr
call genlofr
call getufr
call genufrp
call genwfnr(151,tq0bz)
if (mpi_grid_root()) then
  fname="wfnrkp.hdf5"
  call hdf5_create_file(fname)
  call hdf5_create_group(fname,"/","parameters")
  call hdf5_create_group(fname,"/","kpoints")
endif

if (mpi_grid_side(dims=(/dim_k/))) then
  do i=0,mpi_grid_size(dim_k)-1
    if (mpi_grid_x(dim_k).eq.i) then
      do ikloc=1,nkptnrloc
        ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
        write(kname,'(I10.10)')ik
        call hdf5_create_group(fname,"/kpoints",kname)
        call hdf5_write(fname,"/kpoints/"//trim(kname),"wfsvmt",&
          wfsvmtnrloc(1,1,1,1,1,ikloc),(/lmmaxvr,nufrmax,natmtot,nspinor,nstsv/))
        call hdf5_write(fname,"/kpoints/"//trim(kname),"wfsvit",&
          wfsvitnrloc(1,1,1,ikloc),(/ngkmax,nspinor,nstsv/))
        call hdf5_write(fname,"/kpoints/"//trim(kname),"evalsv",&
          evalsvnr(1,ik),(/nstsv/))
        call hdf5_write(fname,"/kpoints/"//trim(kname),"occsv",&
          occsvnr(1,ik),(/nstsv/))
        if (allocated(pmatnrloc)) then
          call hdf5_write(fname,"/kpoints/"//trim(kname),"pmat",&
            pmatnrloc(1,1,1,ikloc),(/3,nstsv,nstsv/))
        endif
      enddo
    endif
    call mpi_grid_barrier(dims=(/dim_k/))
  enddo
endif
return
end