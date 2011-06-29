subroutine sic_update_umtrx
use modmain
use mod_sic
implicit none
integer n,m,ikloc,i,j,n1,n2,j1,j2,vl(3),ik
real(8) tot_diff
complex(8), allocatable :: um(:,:),um1(:,:),um2(:,:)
allocate(um(sic_wantran%nwan,sic_wantran%nwan))
allocate(um1(nwantot,nwantot))
allocate(um2(nwantot,nwantot))
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  tot_diff=0.d0
  call sic_exp_grad_u(vkcnr(:,ik),sic_umtrx_eps,tot_diff,um)
  if (debug_level.ge.4) then
    call dbg_open_file
    write(fdbgout,'("k-point : ",I4)')ik
    do n1=1,nwantot
      write(fdbgout,'(255F12.6)')(abs(um(n1,n2)),n2=1,nwantot)
    enddo
    call dbg_close_file
  endif
  um2=zzero
  do i=1,nwantot
    um2(i,i)=zone
  enddo
  do j1=1,sic_wantran%nwan
    n1=sic_wantran%iwan(j1)
    do j2=1,sic_wantran%nwan
      n2=sic_wantran%iwan(j2)
      um2(n1,n2)=um(j1,j2)
    enddo
  enddo
  um1=zzero
  do n1=1,nwantot
    do n2=1,nwantot
      do n=1,nwantot
        um1(n1,n2)=um1(n1,n2)+sic_wan_umtrx(n1,n,ikloc)*um2(n,n2)
      enddo
    enddo
  enddo
  sic_wan_umtrx(:,:,ikloc)=um1(:,:)
enddo
deallocate(um,um1,um2)
return
end
