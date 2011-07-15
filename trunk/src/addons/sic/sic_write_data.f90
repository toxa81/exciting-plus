subroutine sic_write_data
use modmain
use mod_sic
use mod_hdf5
use mod_linresp
use mod_wannier
implicit none
integer n,ispn,j,i,ikloc,ik
character*12 c1,c2
character*100 path
!
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
    enddo
  enddo
  call hdf5_create_group("sic.hdf5","/","kpoint")
  call hdf5_write("sic.hdf5","/kpoint","nkpt",nkpt)
  do ik=1,nkpt
    path="/kpoint"
    write(c1,'(I4.4)')ik
    call hdf5_create_group("sic.hdf5",path,trim(adjustl(c1)))
    do n=1,nwantot
      path="/kpoint/"//trim(adjustl(c1))
      write(c2,'("n",I4.4)')n
      call hdf5_create_group("sic.hdf5",path,trim(adjustl(c2)))   
    enddo
  enddo
  call hdf5_write("sic.hdf5","/","nwt",sic_wantran%nwt)
  call hdf5_write("sic.hdf5","/","iwt",sic_wantran%iwt(1,1),&
    (/5,sic_wantran%nwt/))
  call hdf5_write("sic.hdf5","/","vme",sic_vme(1),(/sic_wantran%nwt/))
  call hdf5_write("sic.hdf5","/","e0",sic_wan_e0(1),(/nwantot/))
  call hdf5_write("sic.hdf5","/","sic_energy_tot",sic_energy_tot)
  call hdf5_write("sic.hdf5","/","sic_energy_pot",sic_energy_pot)
  call hdf5_write("sic.hdf5","/","sic_energy_kin",sic_energy_kin)  
  do n=1,nwantot
    j=sic_wantran%idxiwan(n)
    if (j.gt.0) then
      do ispn=1,nspinor
        write(c1,'("n",I4.4)')n
        write(c2,'("s",I4.4)')ispn
        path="/wann/"//trim(adjustl(c1))//"/"//trim(adjustl(c2))       
        call hdf5_write("sic.hdf5",path,"wlm",s_wlm(1,1,ispn,j),&
          (/lmmaxwan,s_nr/))
        call hdf5_write("sic.hdf5",path,"wvlm",s_wvlm(1,1,ispn,j),&
          (/lmmaxwan,s_nr/))
      enddo
    endif
  enddo
endif
if (mpi_grid_side(dims=(/dim_k/))) then
  do i=0,mpi_grid_dim_size(dim_k)-1
    if (mpi_grid_dim_pos(dim_k).eq.i) then
      do ikloc=1,nkptloc
        ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
        write(c1,'(I4.4)')ik
        do n=1,nwantot
          j=sic_wantran%idxiwan(n)
          if (j.gt.0) then
            write(c2,'("n",I4.4)')n
            path="/kpoint/"//trim(adjustl(c1))//"/"//trim(adjustl(c2))
            call hdf5_write("sic.hdf5",path,"wkmt",s_wkmt(1,1,1,1,j,ikloc),&
              (/nrmtmax,lmmaxapw,natmtot,nspinor/))
            call hdf5_write("sic.hdf5",path,"wkit",s_wkit(1,1,j,ikloc),&
              (/ngkmax,nspinor/))
          endif
        enddo !n
      enddo !ikloc
    end if
    call mpi_grid_barrier(dims=(/dim_k/))
  end do
endif

return
end
