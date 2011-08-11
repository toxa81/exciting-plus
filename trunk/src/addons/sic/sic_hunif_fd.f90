subroutine sic_hunif_fd(ikloc,nmatp,hunif,om)
use modmain
use mod_sic
use mod_wannier
! arguments
implicit none
integer, intent(in) :: ikloc
integer, intent(in) :: nmatp
complex(8), intent(inout) :: hunif(nmatp,nmatp,nspinor)
complex(8), intent(inout) :: om(nmatp,nmatp)
! local variables
integer i,ik,vtrl(3),n1,n2,i1,i2,j1,j2,ispn
real(8) vtrc(3),vn
complex(8) expikt
complex(8), allocatable :: vk(:,:)
character*500 msg,fname
logical, parameter :: tcheckherm=.false.
!
ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
call timer_start(t_sic_hunif)
sic_wan_h0k(:,:,ikloc)=zzero
! compute V_{nn'}(k)
allocate(vk(sic_wantran%nwan,sic_wantran%nwan))
vk=zzero
do i=1,sic_wantran%nwt
  n1=sic_wantran%iwt(1,i)
  j1=sic_wantran%idxiwan(n1)
  n2=sic_wantran%iwt(2,i)
  j2=sic_wantran%idxiwan(n2)
  vtrl(:)=sic_wantran%iwt(3:5,i)
  vtrc(:)=vtrl(1)*avec(:,1)+vtrl(2)*avec(:,2)+vtrl(3)*avec(:,3)
  expikt=exp(zi*dot_product(vkc(:,ik),vtrc(:)))
  vk(j1,j2)=vk(j1,j2)+expikt*sic_vme(i)
enddo
! setup unified Hamiltonian
do ispn=1,nspinor
  do i1=1,nmatp
    do i2=1,nmatp
      do j1=1,sic_wantran%nwan
        n1=sic_wantran%iwan(j1)
        vn=dreal(sic_vme(sic_wantran%iwtidx(n1,n1,0,0,0)))
        do j2=1,sic_wantran%nwan
! -\sum_{\alpha,\alpha'} ( P_{\alpha}V_{\alpha}P_{\alpha'} +
!                           P_{\alpha'}V_{\alpha}P_{\alpha} )
          hunif(i1,i2,ispn)=hunif(i1,i2,ispn)-&
            vk(j1,j2)*dconjg(sic_wb(j1,i1,ispn,ikloc))*&
            sic_wb(j2,i2,ispn,ikloc)-dconjg(vk(j1,j2))*&
            dconjg(sic_wb(j2,i1,ispn,ikloc))*sic_wb(j1,i2,ispn,ikloc)
        enddo !j2
! \sum_{alpha} P_{\alpha} V_{\alpha} P_{\alpha}
        hunif(i1,i2,ispn)=hunif(i1,i2,ispn)+&
          dconjg(sic_wb(j1,i1,ispn,ikloc))*&
          sic_wb(j1,i2,ispn,ikloc)*vn
! \sum_{\alpha} P_{\alpha}V_{\alpha}+V_{\alpha}P_{\alpha}
        hunif(i1,i2,ispn)=hunif(i1,i2,ispn)+&
          dconjg(sic_wb(j1,i1,ispn,ikloc))*sic_wvb(j1,i2,ispn,ikloc)+&
          dconjg(sic_wvb(j1,i1,ispn,ikloc))*sic_wb(j1,i2,ispn,ikloc)
      enddo !j1
    enddo !i2 
  enddo !i1
enddo !ispn
deallocate(vk)
call timer_stop(t_sic_hunif)
return
end subroutine

