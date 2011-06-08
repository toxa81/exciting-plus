subroutine setovl(ngp,nmatp,igpig,apwalm,o)
use modmain
implicit none
integer, intent(in) :: ngp
integer, intent(in) :: nmatp
integer, intent(in) :: igpig(ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(out) :: o(nmatp,nmatp)
complex(8), allocatable :: zm1(:,:)

complex(8) zt1
integer is,ia,ias,ig,naa
integer l1,m1,lm1,io1
integer i,j,ilo1,ilo2
integer iv(3)

call timer_start(t_seceqnfv_setup_o)
call timer_start(t_seceqnfv_setup_o_mt)
allocate(zm1(apwordmax*lmmaxapw*natmtot,ngkmax))
zm1=zzero
naa=0
do ias=1,natmtot
  is=ias2is(ias)
  ia=ias2ia(ias)
  do l1=0,lmaxmat
    do m1=-l1,l1
      lm1=idxlm(l1,m1)
      do io1=1,apword(l1,is)
         naa=naa+1
         zm1(naa,:)=apwalm(:,io1,lm1,ias)
      enddo
    enddo
  enddo
enddo
!----------------------!
!     APW-APW term     !
!----------------------!
call zgemm('C','N',ngp,ngp,naa,zone,zm1,apwordmax*lmmaxapw*natmtot,&
  zm1,apwordmax*lmmaxapw*natmtot,zone,o(1,1),nmatp)
deallocate(zm1)
do ias=1,natmtot
  is=ias2is(ias)
  ia=ias2ia(ias)
!----------------------!
!     APW-APW term     !
!----------------------!
!  do l1=0,lmaxmat
!    do m1=-l1,l1
!      lm1=idxlm(l1,m1)
!      do io1=1,apword(l1,is)
!        do ig=1,ngp
!          call zaxpy(ig,dconjg(apwalm(ig,io1,lm1,ias)),apwalm(:,io1,lm1,ias),1,o(:,ig),1)
!        enddo
!      enddo !io1
!    enddo !m1
!  enddo !l1
!---------------------!
!     APW-lo term     !
!---------------------!  
  do ilo1=1,nlorb(is)
    l1=lorbl(ilo1,is)
    do m1=-l1,l1
      lm1=idxlm(l1,m1)
      i=ngp+idxlo(lm1,ilo1,ias)
      do io1=1,apword(l1,is)
        zt1=zone*oalo(io1,ilo1,ias)
        do ig=1,ngp
          o(ig,i)=o(ig,i)+dconjg(apwalm(ig,io1,lm1,ias))*oalo(io1,ilo1,ias)
        end do
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
call timer_stop(t_seceqnfv_setup_o_mt)
call timer_start(t_seceqnfv_setup_o_it)
!---------------------!
!     interstitial    !
!---------------------!
do j=1,ngp
  do i=1,j
    iv(:)=ivg(:,igpig(i))-ivg(:,igpig(j))
    ig=ivgig(iv(1),iv(2),iv(3))
    o(i,j)=o(i,j)+cfunig(ig)
  end do
end do
call timer_stop(t_seceqnfv_setup_o_it)
call timer_stop(t_seceqnfv_setup_o)
return
end

