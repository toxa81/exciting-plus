subroutine sic_wan_ene
use modmain
use mod_sic
use mod_nrkp
implicit none
integer j,n,ik,ikloc
!
sic_wan_e0=0.d0
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    sic_wan_e0(n)=sic_wan_e0(n)+sic_wan_h0k(j,j,ikloc)*wkpt(ik)
  enddo
enddo
call mpi_grid_reduce(sic_wan_e0(1),nwantot,dims=(/dim_k/),all=.true.)
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  if (mpi_grid_root()) then
    write(60,'(" n : ",I4,"   LDA energy : ",F12.6)')n,sic_wan_e0(n)
  endif
enddo
return
end subroutine
