subroutine wfsvplot
use modmain
implicit none
integer ik,ist,ikloc,j
complex(8), allocatable :: evecfv(:,:)
complex(8), allocatable :: evecsv(:,:)
complex(8), allocatable :: wfsvmt(:,:,:,:,:)
complex(8), allocatable :: wfsvit(:,:,:)
complex(8), allocatable :: apwalm(:,:,:,:)
integer ip,np,ip1,ip2,ip3,nploc,iploc,i
real(8), allocatable :: vrc(:,:)
complex(8), allocatable :: zf(:,:)
real(8) v1(3),v2(3),v3(3),t1,t2,t3
real(8), allocatable :: vgkc1(:,:)

call init0
call init1
call readstate
call linengy
call genapwfr
call genlofr
call geturf
call genurfprod
allocate(evecfv(nmatmax,nstfv))
allocate(evecsv(nstsv,nstsv))
ik=kstlist(1,1)
ist=kstlist(2,1)
ikloc=mpi_grid_map(nkpt,dim_k,glob=ik,x=j)
allocate(wfsvmt(lmmaxvr,nrfmax,natmtot,nspinor,nstsv))
allocate(wfsvit(ngkmax,nspinor,nstsv))
allocate(vgkc1(3,ngkmax))
if (mpi_grid_x(dim_k).eq.j) then
  call getevecfv(vkl(1,ik),vgkl(1,1,1,ikloc),evecfv)
  call getevecsv(vkl(1,ik),evecsv)
  allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
  call match(ngk(1,ik),gkc(:,1,ikloc),tpgkc(:,:,1,ikloc),        &
    sfacgk(:,:,1,ikloc),apwalm)
  call genwfsvmt(lmaxvr,lmmaxvr,ngk(1,ik),evecfv,evecsv,apwalm,wfsvmt)
  call genwfsvit(ngk(1,ik),evecfv,evecsv,wfsvit)
  vgkc1(:,:)=vgkc(:,:,1,ikloc)
endif
call mpi_grid_bcast(wfsvmt(1,1,1,1,1),lmmaxvr*nrfmax*natmtot*nspinor*nstsv,&
  dims=(/dim_k/),root=(/j/))
call mpi_grid_bcast(wfsvit(1,1,1),ngkmax*nspinor*nstsv,dims=(/dim_k/),&
  root=(/j/))
call mpi_grid_bcast(vgkc1(1,1),3*ngkmax,dims=(/dim_k/),root=(/j/))  

allocate(vrc(3,np3d(1)*np3d(2)*np3d(3)))
v1(:)=vclp3d(:,2) 
v2(:)=vclp3d(:,3) 
v3(:)=vclp3d(:,4) 
if (mpi_grid_root()) then
  write(*,*)
  write(*,'("Info(plot3d): cartesian vectors of the plane : ")')
  write(*,'("  v1     : ",3G18.10)')v1
  write(*,'("  v2     : ",3G18.10)')v2
  write(*,'("  v3     : ",3G18.10)')v3
  write(*,'("  center : ",3G18.10)')vclp3d(:,1)
  write(*,*)
endif
ip=0
do ip1=0,np3d(1)-1
  t1=dble(ip1)/dble(np3d(1))
  do ip2=0,np3d(2)-1
    t2=dble(ip2)/dble(np3d(2))
    do ip3=0,np3d(3)-1
      t3=dble(ip3)/dble(np3d(3))
      ip=ip+1
!      vrc(:,ip)=vclp3d(:,1)-0.5d0*(v1(:)+v2(:)+v3(:))+t1*v1(:)+t2*v2(:)+t3*v3(:)
      vrc(:,ip)=vclp3d(:,1)+t1*v1(:)+t2*v2(:)+t3*v3(:)
    end do
  end do
end do
np=ip
allocate(zf(nspinor,np))
zf=zzero
nploc=mpi_grid_map(np,dim_k)
do iploc=1,nploc
  ip=mpi_grid_map(np,dim_k,loc=iploc)
  call wfsv_val(vrc(:,ip),vkc(:,ik),vgkc1,ngk(1,ik),&
    wfsvmt(1,1,1,1,ist),wfsvit(1,1,ist),zf(:,ip))
enddo
call mpi_grid_reduce(zf(1,1),nspinor*np)
if (mpi_grid_root()) then
  open(70,file="WFSV.dx",form="FORMATTED",status="REPLACE")
  write(70,'("object 1 class gridpositions counts",3I4)')np3d(:)
  write(70,'("origin ",3G18.10)')vclp3d(:,1)
  do i=1,3
    write(70,'("delta ",3G18.10)')vclp3d(:,i+1)/np3d(i)
  enddo
  write(70,'("object 2 class gridconnections counts",3I4)')np3d(:)
  write(70,'("object 3 class array type float rank 1 shape 1 items ",&
    &I8," lsb ieee data 4")')np
  write(70,'("object ""regular positions regular connections"" class field")')
  write(70,'("component ""positions"" value 1")')
  write(70,'("component ""connections"" value 2")')
  write(70,'("component ""data"" value 3")')
  write(70,'("end")')
  close(70)
  open(70,file="WFSV.dx",status="OLD",form="UNFORMATTED",position="APPEND")
  write(70)(sngl(sum(abs(zf(:,ip))**2)),ip=1,np)
  close(70)
!  open(160,file="WFSV.OUT",form="FORMATTED",status="REPLACE")
!  write(160,'(3I6," : grid size")') np3d(:)
!  do ip=1,np
!    write(160,'(7G18.10)') vrc(:,ip),sum(abs(zf(:,ip))**2)
!  end do
!  close(160)
endif
return
end