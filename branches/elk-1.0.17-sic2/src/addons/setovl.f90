subroutine setovl(ngp,nmatp,igpig,apwalm,o)
use modmain
implicit none
integer, intent(in) :: ngp
integer, intent(in) :: nmatp
integer, intent(in) :: igpig(ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(out) :: o(nmatp,nmatp)

complex(8) zt1
integer is,ia,ias,ig
integer l1,m1,lm1,io1
integer i,j,ilo1,ilo2
integer natmtotloc,iasloc
integer iv(3)

call timer_start(t_seceqnfv_setup_o)
call timer_start(t_seceqnfv_setup_o_mt)
natmtotloc=mpi_grid_map(natmtot,dim2)
do iasloc=1,natmtotloc
  ias=mpi_grid_map(natmtot,dim2,loc=iasloc)
  is=ias2is(ias)
  ia=ias2ia(ias)
!----------------------!
!     APW-APW term     !
!----------------------!
  do l1=0,lmaxmat
    do m1=-l1,l1
      lm1=idxlm(l1,m1)
      do io1=1,apword(l1,is)
        do ig=1,ngp
          call zaxpy(ig,dconjg(apwalm(ig,io1,lm1,ias)),apwalm(:,io1,lm1,ias),1,o(:,ig),1)
        enddo
      enddo !io1
    enddo !m1
  enddo !l1
!---------------------!
!     APW-lo term     !
!---------------------!  
  do ilo1=1,nlorb(is)
    l1=lorbl(ilo1,is)
    do m1=-l1,l1
      lm1=idxlm(l1,m1)
      i=ngp+idxlo(lm1,ilo1,ias)
      do io1=1,apword(l1,is)
        zt1=dcmplx(oalo(io1,ilo1,ias),0.d0)
        call zaxpy(ngp,zt1,apwalm(:,io1,lm1,ias),1,o(:,i),1)
      end do
    end do
  end do
!--------------------!
!     lo-lo term     !
!--------------------!
  do ilo1=1,nlorb(is)
    l1=lorbl(ilo1,is)
    do ilo2=1,nlorb(is)
      if (lorbl(ilo2,is).eq.l1) then
        do m1=-l1,l1
          lm1=idxlm(l1,m1)
          i=ngp+idxlo(lm1,ilo1,ias)
          j=ngp+idxlo(lm1,ilo2,ias)
          if (i.le.j) then
            o(i,j)=o(i,j)+dcmplx(ololo(ilo1,ilo2,ias),0.d0)
          end if
        end do
      end if
    end do
  end do
enddo
call mpi_grid_reduce(o(1,1),nmatp*nmatp,dims=(/dim2/))
call timer_stop(t_seceqnfv_setup_o_mt)
call timer_start(t_seceqnfv_setup_o_it)
!---------------------!
!     interstitial    !
!---------------------!
do j=1,ngp
  do i=1,j
    iv(:)=ivg(:,igpig(i))-ivg(:,igpig(j))
    ig=ivgig(iv(1),iv(2),iv(3))
    o(i,j)=o(i,j)+dconjg(cfunig(ig))
  end do
end do
call timer_stop(t_seceqnfv_setup_o_it)
o=dconjg(o)
call timer_stop(t_seceqnfv_setup_o)
return
end

