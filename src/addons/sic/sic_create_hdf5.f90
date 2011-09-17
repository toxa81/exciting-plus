subroutine sic_create_hdf5
use modmain
use mod_sic
use mod_hdf5
implicit none
!
integer n,ik,j
character*10 c1
character*100 path
!
if (wproc) then
  call hdf5_create_file("sic.hdf5")
  call hdf5_create_group("sic.hdf5","/","wannier_functions")
  call hdf5_write("sic.hdf5","/","nkpt",nkpt)
  call hdf5_write("sic.hdf5","/","lmmaxwan",lmmaxwan)
  call hdf5_write("sic.hdf5","/","nspinor",nspinor)
  call hdf5_write("sic.hdf5","/","s_nr_min",s_nr_min)
  call hdf5_write("sic.hdf5","/","s_nr",s_nr)
  call hdf5_write("sic.hdf5","/","s_r",s_r(1),(/s_nr/))
  call hdf5_write("sic.hdf5","/","nwt",sic_wantran%nwt)
  call hdf5_write("sic.hdf5","/","iwt",sic_wantran%iwt(1,1),&
    (/5,sic_wantran%nwt/))
  do n=1,nwantot
    j=sic_wantran%idxiwan(n)
    if (j.gt.0) then
      path="/wannier_functions"
      write(c1,'(I4.4)')n
      call hdf5_create_group("sic.hdf5",path,trim(adjustl(c1)))
      path=trim(path)//"/"//trim(adjustl(c1))
      call hdf5_create_group("sic.hdf5",path,"kpoints") 
      path=trim(path)//"/"//"kpoints"
      do ik=1,nkpt
        write(c1,'(I4.4)')ik
        call hdf5_create_group("sic.hdf5",path,trim(adjustl(c1)))
      enddo
    endif
  enddo
endif
return
end subroutine
