subroutine sumchi0wan_k(ikloc,w,chi0wan_k)
use modmain
implicit none
integer, intent(in) :: ikloc
complex(8), intent(in) :: w
complex(8), intent(out) :: chi0wan_k(nmegqwan,nmegqwan)
integer ik,jk,i1,i,ibloc,xloc,ist1,ist2,n,n1,n2
complex(8) wt
complex(8), allocatable :: zv1(:)
allocate(zv1(nmegqwan))
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
jk=idxkq(1,ik)
do i1=1,nmegqblhwan(ikloc)
  i=imegqblhwan(i1,ikloc)
  ibloc=mpi_grid_map(nmegqblh(ikloc),dim_b,glob=i,x=xloc)
  if (xloc.eq.mpi_grid_x(dim_b)) then
    ist1=bmegqblh(1,i,ikloc)
    ist2=bmegqblh(2,i,ikloc)
    do n=1,nmegqwan
      n1=bmegqwan(1,n)
      n2=bmegqwan(2,n)
      zv1(n)=wann_c(n1,ist1,ikloc)*dconjg(wann_c(n2,ist2,ikloc+nkptnrloc))
    enddo
    wt=(lr_occsvnr(ist1,ik)-lr_occsvnr(ist2,jk))/(lr_evalsvnr(ist1,ik) - &
        lr_evalsvnr(ist2,jk)+w)
    call zgerc(nmegqwan,nmegqwan,wt,zv1,1,zv1,1,chi0wan_k,nmegqwan)
  endif
enddo !i1
deallocate(zv1)
return
end