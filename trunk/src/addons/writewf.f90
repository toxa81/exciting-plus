subroutine writewf
use modmain
use mod_nrkp
use mod_hdf5
implicit none
character*100 fname,path,path1
character*20 c1,c2
integer ikloc,ik,i,is,ia,ias,ic

call init0
call init1
if (.not.mpi_grid_in()) return
call readstate
call linengy
call genapwfr
call genlofr
call getufr
call genufrp
wproc=mpi_grid_root()
call genwfnr(6,.true.)

fname="wfnrkp.hdf5"
if (mpi_grid_root()) then
  call hdf5_create_file(fname)
  call hdf5_create_group(fname,"/","parameters")
  call hdf5_create_group(fname,"/","species")
  call hdf5_create_group(fname,"/","kpoints")
  call hdf5_write(fname,"/parameters","avec",avec(1,1),(/3,3/))
  call hdf5_write(fname,"/parameters","ainv",ainv(1,1),(/3,3/))
  call hdf5_write(fname,"/parameters","bvec",bvec(1,1),(/3,3/))
  call hdf5_write(fname,"/parameters","omega",omega)
  call hdf5_write(fname,"/parameters","ngridk",ngridk(1),(/3/))
  call hdf5_write(fname,"/parameters","nkptnr",nkptnr)
  call hdf5_write(fname,"/parameters","lmaxapw",lmaxapw)
  call hdf5_write(fname,"/parameters","natmtot",natmtot)
  call hdf5_write(fname,"/parameters","nspecies",nspecies)
  call hdf5_write(fname,"/parameters","nstsv",nstsv)
  call hdf5_write(fname,"/parameters","nspinor",nspinor)
  call hdf5_write(fname,"/parameters","nufrmax",nufrmax)
  call hdf5_write(fname,"/parameters","ngkmax",ngkmax)
  call hdf5_write(fname,"/parameters","nrmtmax",nrmtmax)
  do is=1,nspecies
    write(c1,'(I3.3)')is
    call hdf5_create_group(fname,"/species",c1)
    path="/species/"//trim(adjustl(c1))
    call hdf5_write(fname,path,"nrmt",nrmt(is))
    call hdf5_write(fname,path,"rmt",rmt(is))
    call hdf5_write(fname,path,"spr",spr(1,is),(/nrmt(is)/))
    call hdf5_write(fname,path,"nufr",nufr(0,is),(/lmaxapw+1/))
    call hdf5_write(fname,path,"natoms",natoms(is))
    call hdf5_create_group(fname,path,"atoms")
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      ic=ias2ic(ias)
      write(c2,'(I3.3)')ia
      call hdf5_create_group(fname,trim(path)//"/atoms/",c2)
      path1=trim(path)//"/atoms/"//trim(adjustl(c2))
      call hdf5_write(fname,path1,"atposc",atposc(1,ia,is),(/3/))
      call hdf5_write(fname,path1,"atposl",atposl(1,ia,is),(/3/))
      call hdf5_write(fname,path1,"ufr",ufr(1,0,1,ic),(/nrmtmax,lmaxapw+1,nufrmax/))
    enddo
  enddo
endif
if (mpi_grid_side(dims=(/dim_k/))) then
  do i=0,mpi_grid_dim_size(dim_k)-1
    if (mpi_grid_dim_pos(dim_k).eq.i) then
      do ikloc=1,nkptnrloc
        ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
        write(c1,'(I4.4)')ik
        call hdf5_create_group(fname,"/kpoints",c1)
        path="/kpoints/"//trim(adjustl(c1))
        call hdf5_write(fname,path,"vkl",vklnr(1,ik),(/3/))
        call hdf5_write(fname,path,"vkc",vkcnr(1,ik),(/3/))
        call hdf5_write(fname,path,"evalsv",evalsvnr(1,ik),(/nstsv/))
        call hdf5_write(fname,path,"occsv",occsvnr(1,ik),(/nstsv/))
        call hdf5_write(fname,path,"wfsvmt",wfsvmtnrloc(1,1,1,1,1,ikloc),&
          (/lmmaxapw,nufrmax,natmtot,nspinor,nstsv/))
        call hdf5_write(fname,path,"wfsvit",wfsvitnrloc(1,1,1,ikloc),&
          (/ngkmax,nspinor,nstsv/))
        call hdf5_write(fname,path,"ngk",ngknr(ik))
        call hdf5_write(fname,path,"vgkl",vgklnr(1,1,ikloc),(/3,ngkmax/))
        call hdf5_write(fname,path,"vgkc",vgkcnr(1,1,ikloc),(/3,ngkmax/))
      enddo
    endif
    call mpi_grid_barrier(dims=(/dim_k/))
  enddo
endif
return
end subroutine

