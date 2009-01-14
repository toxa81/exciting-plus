subroutine setovl(ngp,nmatp,igpig,apwalm,o)
use modmain
implicit none
integer, intent(in) :: ngp
integer, intent(in) :: nmatp
integer, intent(in) :: igpig(ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(out) :: o(nmatp,nmatp)

complex(8) zsum,zt1
real(8) t1
integer i1,i2,is,ia,ias,ig
integer l1,m1,lm1,l2,m2,lm2,l3,m3,lm3,io1,io2
integer i,j,ilo1,ilo2
integer iv(3)

!----------------------!
!     APW-APW term     !
!----------------------!
do ig=1,ngp
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do l1=0,lmaxmat
        do io1=1,apword(l1,is)
          do m1=-l1,l1
            lm1=idxlm(l1,m1)
            call zaxpy(ig,dconjg(apwalm(ig,io1,lm1,ias)),apwalm(:,io1,lm1,ias),1,o(:,ig),1)
          enddo !m1
        enddo !io1
      enddo !l1
    enddo !ia
  enddo !is
enddo !ig

do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
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
		  call zaxpy(ngp,zt1,dconjg(apwalm(:,io1,lm1,ias)),1,o(:,i),1)
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
  end do !ia
end do !is

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

return
end

