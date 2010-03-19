subroutine genmegqwan(ikloc,wann_c_jk)
use modmain
implicit none
integer, intent(in) :: ikloc
complex(8), intent(in) :: wann_c_jk(nwann,nstsv,nkptnrloc)
integer ik,i1,ibloc,xloc,ist1,ist2,itr,i,n1,n2,ig,n,j
real(8) vtrc(3)
complex(8) zt1
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
do i1=1,nmegqblhwan(ikloc)
  i=imegqblhwan(i1,ikloc)
  ibloc=mpi_grid_map(nmegqblh(ikloc),dim_b,glob=i,x=xloc)
  if (xloc.eq.mpi_grid_x(dim_b)) then
    ist1=bmegqblh(1,i,ikloc)
    ist2=bmegqblh(2,i,ikloc)
    do j=1,nmegqwan
      n1=imegqwan(1,j)
      n2=imegqwan(2,j)
      vtrc(:)=avec(:,1)*imegqwan(3,j)+&
              avec(:,2)*imegqwan(4,j)+&
              avec(:,3)*imegqwan(5,j)
      zt1=exp(dcmplx(0.d0,-dot_product(vkcnr(:,ik)+vq0rc(:),vtrc(:))))
      do ig=1,ngvecme
        megqwan(j,ig)=megqwan(j,ig)+dconjg(wann_c(n1,ist1,ikloc))*&
          wann_c_jk(n2,ist2,ikloc)*megqblh(ibloc,ig,ikloc)*zt1
      enddo !ig      
    enddo !j
  endif !if (xloc.eq.mpi_grid_x(dim_b))
enddo !i1
return
end