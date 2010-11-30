subroutine genchi0wan_k(ikloc,w,chi0wan_k)
use modmain
use mod_nrkp
use mod_linresp
use mod_wannier
implicit none
integer, intent(in) :: ikloc
complex(8), intent(in) :: w
complex(8), intent(out) :: chi0wan_k(nmegqwan,nmegqwan)
integer ik,jk,i1,i,ist1,ist2
complex(8), allocatable :: wt(:)
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
jk=idxkq(1,ik)
allocate(wt(nmegqblhwanmax))
! TODO: utilize second dimension
wt=zzero
do i1=1,nmegqblhwan(ikloc)
  i=imegqblhwan(i1,ikloc)
  ist1=bmegqblh(1,i,ikloc)
  ist2=bmegqblh(2,i,ikloc)
  wt(i1)=(occsvnr(ist1,ik)-occsvnr(ist2,jk))/(evalsvnr(ist1,ik) - &
      evalsvnr(ist2,jk)+w)
enddo
do i1=1,nmegqwan
  wann_cc2(:,i1)=dconjg(wann_cc(:,i1,ikloc))*wt(:)
enddo
call zgemm('T','N',nmegqwan,nmegqwan,nmegqblhwan(ikloc),zone,&
  wann_cc(1,1,ikloc),nmegqblhwanmax,wann_cc2(1,1),nmegqblhwanmax,&
  zone,chi0wan_k(1,1),nmegqwan)
deallocate(wt)
return
end
