subroutine genchi0blh(ikloc,w,chi0w)
use modmain
use mod_nrkp
implicit none
! arguments
integer, intent(in) :: ikloc
complex(8), intent(in) :: w
complex(8), intent(out) :: chi0w(ngvecme,ngvecme)
! local variables
logical l1
integer i,ist1,ist2,offs,ik,jk,ig
complex(8), allocatable :: wt(:)
logical, external :: bndint
! 
ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
jk=idxkq(1,ik)
offs=nmegqblhloc(2,ikloc)
allocate(wt(nmegqblhlocmax))
wt(:)=zzero
do i=1,nmegqblhloc(1,ikloc)
  ist1=bmegqblh(1,i+offs,ikloc)
  ist2=bmegqblh(2,i+offs,ikloc)
! default : include all interband transitions         
  l1=.true.
! cRPA case : don't include bands in energy window [crpa_e1,crpa_e2]
  if (bndint(ist1,evalsvnr(ist1,ik),chi0_exclude_bands(1),&
      chi0_exclude_bands(2)).and.bndint(ist2,evalsvnr(ist2,jk),&
      chi0_exclude_bands(1),chi0_exclude_bands(2))) l1=.false.
  if (l1) then
    if (abs(occsvnr(ist1,ik)-occsvnr(ist2,jk)).gt.1d-10) then
      wt(i)=(occsvnr(ist1,ik)-occsvnr(ist2,jk))/(evalsvnr(ist1,ik) - &
        evalsvnr(ist2,jk)+w)
    endif
  endif
enddo !i
call papi_timer_start(pt_megqblh2)
do ig=1,ngvecme
  megqblh2(:,ig)=dconjg(megqblh(:,ig,ikloc))*wt(:)
enddo
call papi_timer_stop(pt_megqblh2)
call papi_timer_start(pt_chi0_zgemm)
call zgemm('T','N',ngvecme,ngvecme,nmegqblhloc(1,ikloc),zone,&
  megqblh(1,1,ikloc),nmegqblhlocmax,megqblh2(1,1),nmegqblhlocmax,&
  zone,chi0w(1,1),ngvecme)
call papi_timer_stop(pt_chi0_zgemm)
deallocate(wt)
return
end
