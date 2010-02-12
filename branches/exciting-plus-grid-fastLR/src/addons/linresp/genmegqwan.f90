subroutine genmegqwan(ikloc)
use modmain
implicit none
integer, intent(in) :: ikloc
integer ik,i1,ibloc,xloc,ist1,ist2,itr,i,n1,n2,ig,n
real(8) vtrc(3)
complex(8) zt1
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
do i1=1,nmegqblhwan(ikloc)
  i=imegqblhwan(i1,ikloc)
  ibloc=mpi_grid_map(nmegqblh(ikloc),dim_b,glob=i,x=xloc)
  if (xloc.eq.mpi_grid_x(dim_b)) then
    ist1=bmegqblh(1,i,ikloc)
    ist2=bmegqblh(2,i,ikloc)
    do itr=1,ntrmegqwan
      vtrc(:)=avec(:,1)*itrmegqwan(1,itr)+&
              avec(:,2)*itrmegqwan(2,itr)+&
              avec(:,3)*itrmegqwan(3,itr)
      zt1=exp(dcmplx(0.d0,-dot_product(vkcnr(:,ik)+vq0rc(:),vtrc(:))))
      do ig=1,ngvecme
        do n=1,nmegqwan
          n1=bmegqwan(1,n)
          n2=bmegqwan(2,n)
          megqwan(n,itr,ig)=megqwan(n,itr,ig)+dconjg(wann_c(n1,ist1,ikloc))*&
              wann_c(n2,ist2,ikloc+nkptnrloc)*megqblh(ibloc,ig,ikloc)*zt1
        enddo !n
      enddo !ig      
    enddo !itr
  endif !if (xloc.eq.mpi_grid_x(dim_b))
enddo !i1
return
end