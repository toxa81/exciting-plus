subroutine genmegqwan(iq)
use modmain
use mod_nrkp
use mod_addons_q
use mod_expigqr
use mod_wannier
implicit none
integer, intent(in) :: iq
integer ik,i1,ibloc,xloc,ist1,ist2,i,n1,n2,ig,j,ikloc
complex(8) zt2
complex(8), allocatable :: zm1(:,:)
complex(8), allocatable :: expkqt(:)
real(8), allocatable :: vtrc(:,:)
allocate(vtrc(3,megqwantran%nwt))
do j=1,megqwantran%nwt
  vtrc(:,j)=avec(:,1)*megqwantran%iwt(3,j)+avec(:,2)*megqwantran%iwt(4,j)+&
    &avec(:,3)*megqwantran%iwt(5,j)
enddo
allocate(zm1(nwantot,nwantot))
allocate(expkqt(megqwantran%nwt))
do ikloc=1,nkptnrloc
  ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
  do j=1,megqwantran%nwt
    expkqt(j)=exp(-zi*dot_product(vkcnr(:,ik)+vqc(:,iq),vtrc(:,j)))
  enddo
  do i1=1,nmegqblhwan(ikloc)
    i=imegqblhwan(i1,ikloc)
    ist1=bmegqblh(1,i,ikloc)
    ist2=bmegqblh(2,i,ikloc)
    do n1=1,nwantot
      do n2=1,nwantot
        zm1(n1,n2)=dconjg(wanncnrloc(n1,ist1,ikloc))*wann_c_jk(n2,ist2,ikloc)
      enddo
    enddo
    do j=1,megqwantran%nwt
      n1=megqwantran%iwt(1,j)
      n2=megqwantran%iwt(2,j)
      zt2=zm1(n1,n2)*expkqt(j)
      do ig=1,ngq(iq)
        megqwan(j,ig)=megqwan(j,ig)+zt2*megqblh(i,ig,ikloc)
      enddo !ig 
    enddo !j
  enddo !i1
enddo !ikloc
deallocate(zm1)
deallocate(expkqt)
deallocate(vtrc)
return
end
