subroutine sic_ekin
use modmain
use mod_sic
implicit none
integer ikloc,ik,j,n,i,ist
sic_energy_kin=0.d0
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
  do ist=1,nstsv
    do j=1,sic_wantran%nwan
      n=sic_wantran%iwan(j)
      i=sic_wantran%iwtidx(n,n,0,0,0)
      sic_energy_kin=sic_energy_kin-wkpt(ik)*occsv(ist,ik)*&
        dreal(dconjg(wann_c(n,ist,ikloc))*wann_c(n,ist,ikloc)*sic_vme(i))
    enddo
  enddo !ist
enddo !ikloc
call mpi_grid_reduce(sic_energy_kin,dims=(/dim_k/))
return
end subroutine
