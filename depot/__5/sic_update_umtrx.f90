subroutine sic_update_umtrx
use modmain
use mod_sic
implicit none
integer n,m,ikloc,i,j,n1,n2,j1,j2,vl(3),ik
real(8) vtrc(3),eps
complex(8) z1
real(8), allocatable :: eval(:)
complex(8), allocatable :: gm(:,:),um(:,:),um1(:,:)
!
eps=-0.4d0

allocate(gm(sic_wantran%nwan,sic_wantran%nwan))
allocate(eval(sic_wantran%nwan))
allocate(um(sic_wantran%nwan,sic_wantran%nwan))
allocate(um1(nwantot,nwantot))
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  gm=zzero
  do i=1,sic_wantran%nwt
    n1=sic_wantran%iwt(1,i)
    j1=sic_wantran%idxiwan(n1)
    n2=sic_wantran%iwt(2,i)
    j2=sic_wantran%idxiwan(n2)
    vl(:)=sic_wantran%iwt(3:5,i)
    j=sic_wantran%iwtidx(n2,n1,-vl(1),-vl(2),-vl(3))
    vtrc(:)=vl(1)*avec(:,1)+vl(2)*avec(:,2)+vl(3)*avec(:,3)
    z1=exp(zi*dot_product(vkcnr(:,ik),vtrc(:)))
    gm(j1,j2)=gm(j1,j2)+z1*(dconjg(sic_vme(j))-sic_vme(i))
  enddo
  gm=gm*zi
  do j1=1,sic_wantran%nwan
    do j2=1,sic_wantran%nwan
      if (abs(gm(j1,j2)-dconjg(gm(j2,j1))).gt.1d-10) then
        write(*,*)"error! gm is not hermitian"
        write(*,*)j1,j2
        write(*,*)gm(j1,j2)
        write(*,*)gm(j2,j1)
        call pstop
      endif
    enddo
  enddo
  call diagzhe(sic_wantran%nwan,gm,eval)
  um=zzero
  do j1=1,sic_wantran%nwan
    do j2=1,sic_wantran%nwan
      do j=1,sic_wantran%nwan
        um(j1,j2)=um(j1,j2)+gm(j1,j)*dconjg(gm(j2,j))*exp(-dcmplx(0.d0,eps*eval(j)))
      enddo
    enddo
  enddo

  !do j1=1,sic_wantran%nwan
  !  do j2=1,sic_wantran%nwan
  !    if (abs(um(j1,j2)-dconjg(um(j2,j1))).gt.1d-10) then
  !      write(*,*)"error! um is not hermitian"
  !      write(*,*)j1,j2
  !      write(*,*)um(j1,j2)
  !      write(*,*)um(j2,j1)
  !      call pstop
  !    endif
  !  enddo
  !enddo

  um1=zzero
  do j1=1,sic_wantran%nwan
    n1=sic_wantran%iwan(j1)
    do j2=1,sic_wantran%nwan
      n2=sic_wantran%iwan(j2)
      do j=1,sic_wantran%nwan
        n=sic_wantran%iwan(j)
        um1(n1,n2)=um1(n1,n2)+sic_wan_umtrx(n1,n,ikloc)*um(j,j2)
      enddo
    enddo
  enddo
  sic_wan_umtrx(:,:,ikloc)=um1(:,:)
enddo
deallocate(gm,um,um1,eval)
return
end
