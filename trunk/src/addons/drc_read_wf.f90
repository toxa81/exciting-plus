subroutine drc_read_wf
use modmain
use mod_nrkp
use mod_hdf5
use mod_addons_q
implicit none
integer flg_wan,flg_pmat,j,ik,ikloc
character*20 kname
character*100 fname

#ifndef _HDF5_
write(*,'("Error(drc_read_wf): must be compiled with HDF5 support")')
call pstop
#endif

fname="wfnrkp.hdf5"

call gengknr
call hdf5_read(fname,"/parameters","wannier",flg_wan)
call hdf5_read(fname,"/parameters","pmat",flg_pmat)

if (allocated(evalsvnr)) deallocate(evalsvnr)
allocate(evalsvnr(nstsv,nkptnr))
evalsvnr=0.d0
if (allocated(occsvnr)) deallocate(occsvnr)
allocate(occsvnr(nstsv,nkptnr))
occsvnr=0.d0
if (spinpol) then
  if (allocated(spinor_ud)) deallocate(spinor_ud)
  allocate(spinor_ud(2,nstsv,nkptnr))
  spinor_ud=0
endif
if (allocated(wfsvmtnrloc)) deallocate(wfsvmtnrloc)
allocate(wfsvmtnrloc(lmmaxapw,nufrmax,natmtot,nspinor,nstsv,nkptnrloc))
if (allocated(wfsvitnrloc)) deallocate(wfsvitnrloc)
allocate(wfsvitnrloc(ngkmax,nspinor,nstsv,nkptnrloc))
if (flg_wan.eq.1) then
  if (allocated(wanncnrloc)) deallocate(wanncnrloc)
  allocate(wanncnrloc(nwantot,nstsv,nkptnrloc))
endif
if (flg_pmat.eq.1) then
  if (allocated(pmatnrloc)) deallocate(pmatnrloc)
  allocate(pmatnrloc(3,nstsv,nstsv,nkptnrloc))
endif
if (mpi_grid_side(dims=(/dim_k/))) then
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    write(kname,'(I10.10)')ik
    call hdf5_read(fname,"/kpoints/"//trim(kname),"wfsvmt",&
      wfsvmtnrloc(1,1,1,1,1,ikloc),(/lmmaxapw,nufrmax,natmtot,nspinor,nstsv/))
    call hdf5_read(fname,"/kpoints/"//trim(kname),"wfsvit",&
      wfsvitnrloc(1,1,1,ikloc),(/ngkmax,nspinor,nstsv/))
    call hdf5_read(fname,"/kpoints/"//trim(kname),"evalsv",&
      evalsvnr(1,ik),(/nstsv/))
    call hdf5_read(fname,"/kpoints/"//trim(kname),"occsv",&
      occsvnr(1,ik),(/nstsv/))
    if (flg_pmat.eq.1) then
      call hdf5_read(fname,"/kpoints/"//trim(kname),"pmat",&
        pmatnrloc(1,1,1,ikloc),(/3,nstsv,nstsv/))
    endif
    if (flg_wan.eq.1) then
      call hdf5_read(fname,"/kpoints/"//trim(kname),"wannc",&
        wanncnrloc(1,1,ikloc),(/nwantot,nstsv/))          
    endif
    if (spinpol) then
      call hdf5_read(fname,"/kpoints/"//trim(kname),"spinor_ud",&
        spinor_ud(1,1,ik),(/2,nstsv/))
    endif  
  enddo
endif
do ikloc=1,nkptnrloc
  do j=1,nstsv
    call mpi_grid_bcast(wfsvmtnrloc(1,1,1,1,j,ikloc),lmmaxapw*nufrmax*natmtot*nspinor,&
      dims=ortdims((/dim_k/)))
    call mpi_grid_bcast(wfsvitnrloc(1,1,j,ikloc),ngkmax*nspinor,&
      dims=ortdims((/dim_k/)))
  enddo
  if (flg_pmat.eq.1) then
    call mpi_grid_bcast(pmatnrloc(1,1,1,ikloc),3*nstsv*nstsv,&
      dims=ortdims((/dim_k/)))
  endif
  if (flg_wan.eq.1) then
    call mpi_grid_bcast(wanncnrloc(1,1,ikloc),nwantot*nstsv,&
      dims=ortdims((/dim_k/)))
  endif
enddo
call mpi_grid_reduce(evalsvnr(1,1),nstsv*nkptnr,dims=(/dim_k/),&
  all=.true.,side=.true.)
call mpi_grid_bcast(evalsvnr(1,1),nstsv*nkptnr,dims=ortdims((/dim_k/)))
call mpi_grid_reduce(occsvnr(1,1),nstsv*nkptnr,dims=(/dim_k/),&
  all=.true.,side=.true.)
call mpi_grid_bcast(occsvnr(1,1),nstsv*nkptnr,dims=ortdims((/dim_k/)))
if (spinpol) then
  call mpi_grid_reduce(spinor_ud(1,1,1),2*nstsv*nkptnr,dims=(/dim_k/),&
    all=.true.,side=.true.)
  call mpi_grid_bcast(spinor_ud(1,1,1),2*nstsv*nkptnr,dims=ortdims((/dim_k/)))
endif
return
end
