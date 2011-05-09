subroutine sic_wan_ene
use modmain
use mod_sic
use mod_nrkp
implicit none
integer i,j,n,ik,ikloc


!! write eigenvalues/vectors and occupancies to file
!if (mpi_grid_side(dims=(/dim_k/))) then
!  do i=0,mpi_grid_dim_size(dim_k)-1
!    if (mpi_grid_dim_pos(dim_k).eq.i) then
!      do ikloc=1,nkptloc
!        ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
!        !call putevalfv(ik,evalfv(1,1,ikloc))
!        call putevalsv(ik,evalsv(1,ik))
!        call putoccsv(ik,occsv(1,ik))
!        call putevecfv(ik,evecfvloc(1,1,1,ikloc))
!        call putevecsv(ik,evecsvloc(1,1,ikloc))
!      end do
!    end if
!    call mpi_grid_barrier(dims=(/dim_k/))
!  end do
!endif
!
!call genwfnr(-1,.false.)
!! get LDA energies
!sic_wann_e0=0.d0
!do j=1,sic_wantran%nwan
!  n=sic_wantran%iwan(j)
!  i=sic_wantran%iwtidx(n,n,0,0,0)
!  sic_wann_e0(n)=wann_ene(n)-dreal(vwanme(i))
!  if (mpi_grid_root()) then
!    write(60,'(" n : ",I4,"   LDA energy : ",F12.6)')n,sic_wann_e0(n)
!  endif
!enddo

sic_wann_e0=0.d0
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    sic_wann_e0(n)=sic_wann_e0(n)+sic_wann_h0k(j,j,ikloc)*wkpt(ik)
  enddo
enddo
call mpi_grid_reduce(sic_wann_e0(1),nwantot,dims=(/dim_k/),all=.true.)
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  if (mpi_grid_root()) then
    write(60,'(" n : ",I4,"   LDA energy : ",F12.6)')n,sic_wann_e0(n)
  endif
enddo

return
end subroutine
