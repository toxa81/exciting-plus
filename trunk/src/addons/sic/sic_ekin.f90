subroutine sic_ekin
use modmain
use mod_sic
implicit none
integer ikloc,ik,j,n,i,ist,n1,n2,j1,j2,vtrl(3)
real(8) vtrc(3)
complex(8), allocatable :: vk(:,:)
complex(8) expikt

allocate(vk(sic_wantran%nwan,sic_wantran%nwan))

sic_energy_kin=0.d0
do ikloc=1,nkptloc
  ik=mpi_grid_map(nkpt,dim_k,loc=ikloc)
  vk=zzero
  do i=1,sic_wantran%nwt
    n1=sic_wantran%iwt(1,i)
    j1=sic_wantran%idxiwan(n1)
    n2=sic_wantran%iwt(2,i)
    j2=sic_wantran%idxiwan(n2)
    vtrl(:)=sic_wantran%iwt(3:5,i)
    vtrc(:)=vtrl(1)*avec(:,1)+vtrl(2)*avec(:,2)+vtrl(3)*avec(:,3)
    expikt=exp(-zi*dot_product(vkc(:,ik),vtrc(:)))
    vk(j1,j2)=vk(j1,j2)+expikt*dconjg(sic_vme(i))
  enddo
  do ist=1,nstsv
    do j=1,sic_wantran%nwan
      n=sic_wantran%iwan(j)
      !i=sic_wantran%iwtidx(n,n,0,0,0)
      sic_energy_kin=sic_energy_kin-wkpt(ik)*occsv(ist,ik)*&
        &dreal(dconjg(wann_c(n,ist,ikloc))*wann_c(n,ist,ikloc)*vk(j,j))
    enddo
  enddo !ist
enddo !ikloc
deallocate(vk)
call mpi_grid_reduce(sic_energy_kin,dims=(/dim_k/))
return
end subroutine
