subroutine genchi0wan_k(ikloc,w,chi0wan_k)
use modmain
use mod_nrkp
use mod_expigqr
use mod_linresp
implicit none
integer, intent(in) :: ikloc
complex(8), intent(in) :: w
complex(8), intent(out) :: chi0wan_k(megqwantran%nwt,megqwantran%nwt)
integer ik,jk,i1,i,ist1,ist2
complex(8), allocatable :: wt(:)
real(8) t1,t2
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
jk=idxkq(1,ik)
allocate(wt(nmegqblhwanmax))
! TODO: utilize second dimension
wt=zzero
do i1=1,nmegqblhwan(ikloc)
  i=imegqblhwan(i1,ikloc)
  ist1=bmegqblh(1,i,ikloc)
  ist2=bmegqblh(2,i,ikloc)
  t1=occsvnr(ist1,ik)-occsvnr(ist2,jk)
  t2=sign(scissor,t1)
  wt(i1)=t1/(evalsvnr(ist1,ik)-evalsvnr(ist2,jk)-t2+w)
enddo
do i1=1,megqwantran%nwt
  wann_cc2(:,i1)=dconjg(wann_cc(:,i1,ikloc))*wt(:)
enddo
call zgemm('T','N',megqwantran%nwt,megqwantran%nwt,nmegqblhwan(ikloc),zone,&
  wann_cc(1,1,ikloc),nmegqblhwanmax,wann_cc2(1,1),nmegqblhwanmax,&
  zone,chi0wan_k(1,1),megqwantran%nwt)
deallocate(wt)
return
end
