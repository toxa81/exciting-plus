subroutine genmegqwan(iq)
use modmain
use mod_nrkp
use mod_addons_q
implicit none
integer, intent(in) :: iq
integer ik,i1,ibloc,xloc,ist1,ist2,i,n1,n2,ig,j,ikloc,n
complex(8) zt2
complex(8), allocatable :: zm1(:,:)
complex(8), allocatable :: expkqt(:)
real(8), allocatable :: vtrc(:,:)

if (.not.diag_wan_ibt) then
  allocate(vtrc(3,nmegqwan))
  do j=1,nmegqwan
    vtrc(:,j)=avec(:,1)*imegqwan(3,j)+avec(:,2)*imegqwan(4,j)+&
      avec(:,3)*imegqwan(5,j)
  enddo
  allocate(zm1(nwantot,nwantot))
  allocate(expkqt(nmegqwan))
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do j=1,nmegqwan
      expkqt(j)=exp(dcmplx(0.d0,-dot_product(vkcnr(:,ik)+vqc(:,iq),vtrc(:,j))))
    enddo
    do i1=1,nmegqblhwan(ikloc)
      i=imegqblhwan(i1,ikloc)
      ibloc=mpi_grid_map(nmegqblh(ikloc),dim_b,glob=i,x=xloc)
      if (xloc.eq.mpi_grid_x(dim_b)) then
        ist1=bmegqblh(1,i,ikloc)
        ist2=bmegqblh(2,i,ikloc)
        do n1=1,nwantot
          do n2=1,nwantot
            zm1(n1,n2)=dconjg(wanncnrloc(n1,ist1,ikloc))*wann_c_jk(n2,ist2,ikloc)
          enddo
        enddo
        do j=1,nmegqwan
          n1=imegqwan(1,j)
          n2=imegqwan(2,j)
          zt2=zm1(n1,n2)*expkqt(j)
          do ig=1,ngvecme
            megqwan(j,ig)=megqwan(j,ig)+zt2*megqblh(ibloc,ig,ikloc)
          enddo !ig 
        enddo !j
      endif !if (xloc.eq.mpi_grid_x(dim_b))
    enddo !i1
  enddo !ikloc
  deallocate(zm1)
  deallocate(expkqt)
  deallocate(vtrc)
else
  allocate(zm1(nwantot,1))
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    do i1=1,nmegqblhwan(ikloc)
      i=imegqblhwan(i1,ikloc)
      ibloc=mpi_grid_map(nmegqblh(ikloc),dim_b,glob=i,x=xloc)
      if (xloc.eq.mpi_grid_x(dim_b)) then
        ist1=bmegqblh(1,i,ikloc)
        ist2=bmegqblh(2,i,ikloc)
        do n=1,nwantot
          zm1(n,1)=dconjg(wanncnrloc(n,ist1,ikloc))*wann_c_jk(n,ist2,ikloc)
        enddo
        do n=1,nwantot
          j=idxmegqwan(n,n,0,0,0)
          do ig=1,ngvecme
            megqwan(j,ig)=megqwan(j,ig)+zm1(n,1)*megqblh(ibloc,ig,ikloc)
          enddo !ig 
        enddo !j
      endif !if (xloc.eq.mpi_grid_x(dim_b))
    enddo !i1
  enddo !ikloc
  deallocate(zm1)
endif
return
end
