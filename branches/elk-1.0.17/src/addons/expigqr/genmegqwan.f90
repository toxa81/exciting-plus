subroutine genmegqwan(ikloc,expkqt)
use modmain
use mod_nrkp
implicit none
integer, intent(in) :: ikloc
complex(8), intent(in) :: expkqt(nmegqwan,nkptnr)
integer ik,i1,ibloc,xloc,ist1,ist2,i,n1,n2,ig,j
complex(8) zt2
complex(8), allocatable :: zm1(:,:)
allocate(zm1(nwann,nwann))
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
do i1=1,nmegqblhwan(ikloc)
  i=imegqblhwan(i1,ikloc)
  ibloc=mpi_grid_map(nmegqblh(ikloc),dim_b,glob=i,x=xloc)
  if (xloc.eq.mpi_grid_x(dim_b)) then
    ist1=bmegqblh(1,i,ikloc)
    ist2=bmegqblh(2,i,ikloc)
    zm1=zzero
    do n1=1,nwann
      do n2=1,nwann
        zm1(n1,n2)=dconjg(wanncnrloc(n1,ist1,ikloc))*wann_c_jk(n2,ist2,ikloc)
      enddo
    enddo
    do j=1,nmegqwan
      n1=imegqwan(1,j)
      n2=imegqwan(2,j)
      zt2=zm1(n1,n2)*expkqt(j,ik)
      do ig=1,ngvecme
        megqwan(j,ig)=megqwan(j,ig)+zt2*megqblh(ibloc,ig,ikloc)
      enddo !ig      
    enddo !j
  endif !if (xloc.eq.mpi_grid_x(dim_b))
enddo !i1
deallocate(zm1)
return
end