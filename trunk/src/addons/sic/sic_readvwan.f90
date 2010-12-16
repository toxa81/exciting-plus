subroutine sic_readvwan
use modmain
use mod_sic
use mod_hdf5
use mod_linresp
use mod_wannier
implicit none
integer it,n,ispn,nwt,j
character*20 c1,c2,c3
character*100 path
logical exist
complex(8), allocatable :: fmt(:,:,:)
complex(8), allocatable :: fir(:)

! retutn if wannier functions and potential are initialized
if (tsic_wv) return
inquire(file="sic.hdf5",exist=exist)
if (.not.exist) return

call hdf5_read("sic.hdf5","/","nwt",nwt)
if (nwt.ne.sic_wantran%nwt) then
  write(*,'("Error(sic_readvwan) wrong number of Wannier transitions")')
  write(*,'("  nwt read from file : ",I6)')nwt
  write(*,'("  sic_wantran%nwt : ",I6)')sic_wantran%nwt
  call pstop
endif
call hdf5_read("sic.hdf5","/","vwanme",vwanme(1),(/sic_wantran%nwt/))
call hdf5_read("sic.hdf5","/","sic_etot_correction",sic_etot_correction)
call hdf5_read("sic.hdf5","/","sic_epot",sic_epot)

allocate(fmt(lmmaxvr,nrmtmax,natmtot))
allocate(fir(ngrtot))
do n=1,nwantot
  j=sic_wantran%idxiwan(n)
  if (j.gt.0) then
    do ispn=1,nspinor
      do it=1,sic_orbitals%ntr
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
        call sic_copy_mt_z(.true.,lmmaxvr,fmt,sic_orbitals%wvmt(1,1,it,ispn,j))
        call sic_copy_ir_z(.true.,fir,sic_orbitals%wvir(1,it,ispn,j))
        if (mpi_grid_root()) then
          call hdf5_read("sic.hdf5",path,"wanmt",fmt(1,1,1),&
            (/lmmaxvr,nrmtmax,natmtot/))
          call hdf5_read("sic.hdf5",path,"wanir",fir(1),(/ngrtot/))
        endif
        call mpi_grid_bcast(fmt(1,1,1),lmmaxvr*nrmtmax*natmtot)
        call mpi_grid_bcast(fir(1),ngrtot)
        call sic_copy_mt_z(.true.,lmmaxvr,fmt,sic_orbitals%wanmt(1,1,it,ispn,j))
        call sic_copy_ir_z(.true.,fir,sic_orbitals%wanir(1,it,ispn,j))
      enddo
    enddo !ispn
  endif
enddo
deallocate(fmt,fir)
tsic_wv=.true.
return
end
