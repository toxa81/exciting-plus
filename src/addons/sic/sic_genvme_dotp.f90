subroutine sic_genvme_dotp(t0only)
use modmain
use mod_sic
implicit none
logical, intent(in) :: t0only
integer nwtloc,iloc,i,n,j,n1,j1,vl(3)
real(8) pos1(3),pos2(3)
! compute matrix elements of SIC potential
!  sic_vme(i) = <(W*V)_n|W_{n1,T}>
!  i = {n,n1,T}
sic_vme=zzero
if (t0only) then
  nwtloc=mpi_grid_map(sic_wantran%nwt0,dim_k)
else
  nwtloc=mpi_grid_map(sic_wantran%nwt,dim_k)
endif
do iloc=1,nwtloc
  if (t0only) then
    j=mpi_grid_map(sic_wantran%nwt0,dim_k,loc=iloc)  
    i=sic_wantran%iwt0(j)
  else
    i=mpi_grid_map(sic_wantran%nwt,dim_k,loc=iloc)
  endif
  n=sic_wantran%iwt(1,i)
  j=sic_wantran%idxiwan(n)
  n1=sic_wantran%iwt(2,i)
  j1=sic_wantran%idxiwan(n1)
  vl(:)=sic_wantran%iwt(3:5,i)
  pos1(:)=wanpos(:,n)
  pos2(:)=wanpos(:,n1)+vl(1)*avec(:,1)+vl(2)*avec(:,2)+vl(3)*avec(:,3)
  sic_vme(i)=s_spinor_dotp(pos1,pos2,s_wvlm(1,1,1,j),s_wlm(1,1,1,j1))
enddo
call mpi_grid_reduce(sic_vme(1),sic_wantran%nwt,dims=(/dim_k/),all=.true.)
return
end subroutine
