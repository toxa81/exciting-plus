subroutine sumchi0wan_k(ikloc,w,chi0wan_k)
use modmain
implicit none
integer, intent(in) :: ikloc
complex(8), intent(in) :: w
complex(8), intent(out) :: chi0wan_k(nmegqwan,nmegqwan)
integer ik,jk,i1,i,ibloc,xloc,ist1,ist2,n,n1,n2
complex(8), allocatable :: wt(:)
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
jk=idxkq(1,ik)
allocate(wt(nmegqblhwanmax))
wt=zzero
do i1=1,nmegqblhwan(ikloc)
  i=imegqblhwan(i1,ikloc)
  ist1=bmegqblh(1,i,ikloc)
  ist2=bmegqblh(2,i,ikloc)
  wt(i1)=(lr_occsvnr(ist1,ik)-lr_occsvnr(ist2,jk))/(lr_evalsvnr(ist1,ik) - &
      lr_evalsvnr(ist2,jk)+w)
enddo
!do i1=1,nmegqwan
!  wann_cc2(:,i1)=dconjg(wann_cc(:,i1,ikloc))*wt(:)
!enddo
!call zgemm('T','N',nmegqwan,nmegqwan,nmegqblhwan(ikloc),zone,&
!  wann_cc(1,1,ikloc),nmegqblhwanmax,wann_cc2(1,1),nmegqblhwanmax,&
!  zone,chi0wan_k(1,1),nmegqwan)
do i1=1,nmegqblhwan(ikloc)
  call zgerc(nmegqwan,nmegqwan,wt(i1),wann_cc(1,i1,ikloc),1, &
    wann_cc(1,i1,ikloc),1,chi0wan_k,nmegqwan)
enddo !i1
! TODO: utilize second dimension
!do i1=1,nmegqblhwan(ikloc)
!  i=imegqblhwan(i1,ikloc)
!  ist1=bmegqblh(1,i,ikloc)
!  ist2=bmegqblh(2,i,ikloc)
!  do n=1,nmegqwan
!    n1=imegqwan(1,n)
!    n2=imegqwan(2,n)
!    zv1(n)=wann_c(n1,ist1,ikloc)*dconjg(wann_c(n2,ist2,ikloc+nkptnrloc))
!  enddo
!  wt=(lr_occsvnr(ist1,ik)-lr_occsvnr(ist2,jk))/(lr_evalsvnr(ist1,ik) - &
!      lr_evalsvnr(ist2,jk)+w)
!  call zgerc(nmegqwan,nmegqwan,wt,zv1,1,zv1,1,chi0wan_k,nmegqwan)
!enddo !i1
deallocate(wt)
return
end