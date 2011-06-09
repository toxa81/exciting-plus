subroutine sethml(ngp,nmatp,vgpc,igpig,apwalm,h)
use modmain
use mod_timer
implicit none
integer, intent(in) :: ngp
integer, intent(in) :: nmatp
real(8), intent(in) :: vgpc(3,ngkmax)
integer, intent(in) :: igpig(ngkmax)
complex(8), intent(in) :: apwalm(ngkmax,apwordmax,lmmaxapw,natmtot)
complex(8), intent(out) :: h(nmatp,nmatp)

complex(8), allocatable :: zv(:)
complex(8), allocatable :: zm1(:,:)
complex(8), allocatable :: zm2(:,:)
complex(8) zsum,zt1
real(8) t1
integer is,ia,ias,ig,naa
integer l1,m1,lm1,l2,m2,lm2,l3,m3,lm3,io1,io2
integer i,j,ilo1,ilo2
integer iv(3)

call timer_start(t_seceqnfv_setup_h)
call timer_start(t_seceqnfv_setup_h_mt)
allocate(zv(ngp))
allocate(zm1(lmmaxapw*apwordmax*natmtot,ngp))
allocate(zm2(lmmaxapw*apwordmax*natmtot,ngp))
naa=0
do ias=1,natmtot
  is=ias2is(ias)
  ia=ias2ia(ias)
  do l1=0,lmaxapw
    do io1=1,apword(l1,is)
      do m1=-l1,l1
        lm1=idxlm(l1,m1)
        zv=zzero
        do l2=0,lmaxapw
          do m2=-l2,l2
            lm2=idxlm(l2,m2)
            do io2=1,apword(l2,is)
              zsum=0.d0
              do l3=0,lmaxvr
                if (mod(l1+l2+l3,2).eq.0) then
                  do lm3=l3**2+1,(l3+1)**2
                    zsum=zsum+gntyry(lm2,lm3,lm1)*haa(lm3,io2,l2,io1,l1,ias)
                  enddo !m3
                endif !mod(l1+l2+l3,2).eq.0
              enddo !l3
              if (abs(zsum).gt.1.d-14) then
                call zaxpy(ngp,zsum,apwalm(1,io2,lm2,ias),1,zv,1)
              endif
            enddo !io2
          enddo !m2
        enddo !l2
! kinetic surface contribution
        do io2=1,apword(l1,is)
          zt1=(0.5d0*rmt(is)**2)*apwfr(nrmt(is),1,io1,l1,ias)*apwdfr(io2,l1,ias)
          call zaxpy(ngp,zt1,apwalm(1,io2,lm1,ias),1,zv,1)
        enddo !io2
        naa=naa+1
        zm1(naa,:)=zv(:)
        zm2(naa,:)=apwalm(1:ngp,io1,lm1,ias)
      enddo !m1
    enddo !io1
  enddo !l1
enddo
!----------------------!
!     APW-APW term     !
!----------------------!
call zgemm('C','N',ngp,ngp,naa,zone,zm1,lmmaxapw*apwordmax*natmtot,&
  zm2,lmmaxapw*apwordmax*natmtot,zone,h(1,1),nmatp)
deallocate(zm1,zm2,zv)
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
!          call zaxpy(ig,dconjg(zv(ig,lm1,io1,ias)),apwalm(:,io1,lm1,ias),1,h(:,ig),1)
!          call zaxpy(ig,dconjg(apwalm(ig,io1,lm1,ias)),zv(:,lm1,io1,ias),1,h(:,ig),1)
!        enddo !ig
!      enddo !io1
!    enddo !m1
!  enddo !l1
!---------------------!
!     APW-lo term     !
!---------------------!  
  do ilo2=1,nlorb(is)
    l2=lorbl(ilo2,is)
    do m2=-l2,l2
      lm2=idxlm(l2,m2)
      i=ngp+idxlo(lm2,ilo2,ias)
      do l1=0,lmaxapw
        do m1=-l1,l1
          lm1=idxlm(l1,m1)
          do io1=1,apword(l1,is)
            zsum=0.d0
            do l3=0,lmaxvr
              if (mod(l1+l2+l3,2).eq.0) then
                do lm3=l3**2+1,(l3+1)**2
                  zsum=zsum+gntyry(lm1,lm3,lm2)*hloa(lm3,ilo2,io1,l1,ias)
                end do !m3
              end if
            end do !l3
            if (abs(zsum).gt.1.d-14) then
              do ig=1,ngp
                h(ig,i)=h(ig,i)+dconjg(apwalm(ig,io1,lm1,ias))*zsum
              enddo
            endif
          end do !io1
        end do !m1
      end do !l1
    end do !m2
  end do !ilo2
!--------------------!
!     lo-lo term     !
!--------------------!
  do ilo1=1,nlorb(is)
    l1=lorbl(ilo1,is)
    do m1=-l1,l1
      lm1=idxlm(l1,m1)
      i=ngp+idxlo(lm1,ilo1,ias)
      do ilo2=1,nlorb(is)
        l2=lorbl(ilo2,is)
        do m2=-l2,l2
          lm2=idxlm(l2,m2)
          j=ngp+idxlo(lm2,ilo2,ias)
          if (i.le.j) then
            zsum=0.d0
            do l3=0,lmaxvr
              if (mod(l1+l2+l3,2).eq.0) then
                do lm3=l3**2+1,(l3+1)**2
                  zsum=zsum+gntyry(lm1,lm3,lm2)*hlolo(lm3,ilo1,ilo2,ias)
                end do
              end if
            end do
            h(i,j)=h(i,j)+zsum
          end if
        end do
      end do
    end do
  end do
enddo
call timer_stop(t_seceqnfv_setup_h_mt)
!---------------------!
!     interstitial    !
!---------------------!
call timer_start(t_seceqnfv_setup_h_it)
do j=1,ngp
  do i=1,j
    iv(:)=ivg(:,igpig(i))-ivg(:,igpig(j))
    ig=ivgig(iv(1),iv(2),iv(3))
    t1=0.5d0*dot_product(vgpc(:,i),vgpc(:,j))
    h(i,j)=h(i,j)+veffig(ig)+t1*cfunig(ig)
  end do
end do
call timer_stop(t_seceqnfv_setup_h_it)
call timer_stop(t_seceqnfv_setup_h)
return
end

