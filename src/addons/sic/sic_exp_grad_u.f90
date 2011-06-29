subroutine sic_exp_grad_u(vpc,eps,tot_diff,um)
use modmain
use mod_sic
implicit none
real(8), intent(in) :: vpc(3)
real(8), intent(in) :: eps
real(8), intent(out) :: tot_diff
complex(8), intent(out) :: um(sic_wantran%nwan,sic_wantran%nwan)
!
integer i,j,j1,j2,n1,n2,vl(3)
real(8) vtrc(3)
complex(8) expikt
complex(8), allocatable :: gm(:,:)
real(8), allocatable :: eval(:)
!
allocate(gm(sic_wantran%nwan,sic_wantran%nwan))
allocate(eval(sic_wantran%nwan))
gm=zzero
do i=1,sic_wantran%nwt
  n1=sic_wantran%iwt(1,i)
  j1=sic_wantran%idxiwan(n1)
  n2=sic_wantran%iwt(2,i)
  j2=sic_wantran%idxiwan(n2)
  vl(:)=sic_wantran%iwt(3:5,i)
  j=sic_wantran%iwtidx(n2,n1,-vl(1),-vl(2),-vl(3))
  vtrc(:)=vl(1)*avec(:,1)+vl(2)*avec(:,2)+vl(3)*avec(:,3)
  expikt=exp(zi*dot_product(vpc,vtrc(:)))
  gm(j1,j2)=gm(j1,j2)+expikt*(dconjg(sic_vme(j))-sic_vme(i))
enddo
gm=gm*zi
tot_diff=tot_diff+sum(abs(gm))
do j1=1,sic_wantran%nwan
  gm(j1,j1)=gm(j1,j1)+zone
enddo
do j1=1,sic_wantran%nwan
  do j2=1,sic_wantran%nwan
    if (abs(gm(j1,j2)-dconjg(gm(j2,j1))).gt.1d-10) then
      write(*,'("Error(sic_exp_grad_u): gradient matrix is not hermitian")')
      write(*,'("  j1, j2    : ",2I5)')j1,j2
      write(*,'("  gm(j1,j2) : ",2G18.10)')dreal(gm(j1,j2)),dimag(gm(j1,j2))
      write(*,'("  gm(j2,j1) : ",2G18.10)')dreal(gm(j2,j1)),dimag(gm(j2,j1))
      call pstop
    endif
  enddo
enddo
call diagzhe(sic_wantran%nwan,gm,eval)
if (debug_level.ge.4) then
  call dbg_open_file
  write(fdbgout,*)"grad eval=",eval
  call dbg_close_file
endif
um=zzero
! U=e^{eps*G}=e^{-i*i*eps*G}
! H=i*G is a hermitian matrix, so it's decomposition is Z*h*Z^{\dagger}
! finally, e^{-i*eps*H}=Z*e^{-i*eps*h}*Z^{\dagger}
do j1=1,sic_wantran%nwan
  do j2=1,sic_wantran%nwan
    do j=1,sic_wantran%nwan
      um(j1,j2)=um(j1,j2)+gm(j1,j)*dconjg(gm(j2,j))*exp(-zi*eps*(eval(j)-1.d0))
    enddo
  enddo
enddo
deallocate(gm,eval)
return
end subroutine