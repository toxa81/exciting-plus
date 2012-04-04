subroutine sic_write_eval
use modmain
use mod_sic
implicit none
!
integer ik,ikloc,i,n1,n2,j1,j2,vtrl(3)
real(8) vtrc(3)
complex(8) expikt
real(8), allocatable :: e0(:,:)
real(8), allocatable :: e1(:,:)
complex(8), allocatable :: vk(:,:)
!
allocate(e0(sic_wantran%nwan,nkpt))
allocate(e1(sic_wantran%nwan,nkpt))
allocate(vk(sic_wantran%nwan,sic_wantran%nwan))
e0=0.d0
e1=0.d0
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
    expikt=exp(zi*dot_product(vkc(:,ik),vtrc(:)))
    vk(j1,j2)=vk(j1,j2)+expikt*sic_vme(i)
  enddo
  do i=1,sic_wantran%nwan
    e0(i,ik)=dreal(sic_wan_h0k(i,i,ikloc))
    e1(i,ik)=e0(i,ik)+dreal(vk(i,i))
  enddo
enddo
call mpi_grid_reduce(e0(1,1),sic_wantran%nwan*nkpt,dims=(/dim_k/))
call mpi_grid_reduce(e1(1,1),sic_wantran%nwan*nkpt,dims=(/dim_k/))
if (mpi_grid_root()) then
  open(50,file='SIC_EIGVAL'//trim(filext),action='WRITE',form='FORMATTED')
  write(50,'(I6," : nkpt")') nkpt
  do ik=1,nkpt
    write(50,*)
    write(50,'(I6,3G18.10," : k-point, vkl")') ik,vkl(:,ik)
    write(50,'("(wf, e0, e0+v below)")')
    do i=1,sic_wantran%nwan
      write(50,'(I6,2G18.10)') i,e0(i,ik),e1(i,ik)
    enddo
    write(50,*)
  enddo
  close(50)
endif
deallocate(e0,e1,vk)
return
end subroutine
