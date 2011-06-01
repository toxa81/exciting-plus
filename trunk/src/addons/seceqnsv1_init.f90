subroutine seceqnsv1_init
use modmain
use mod_sic
implicit none
real(8), allocatable :: bmt(:,:,:,:)
real(8) cb,t1
integer i,ias,ia,is,ic,l1,l2,l3,m1,m2,m3,lm3,lm1,lm2,io1,io2
integer ir,j1,j2
!
if (.not.spinpol) return

allocate(bmt(lmmaxvr,nrmtmax,natmtot,ndmag))
bmt=0.d0
cb=gfacte/(4.d0*solsc)
do ias=1,natmtot
  is=ias2is(ias)
  ia=ias2ia(ias)
  do i=1,ndmag
    bmt(:,:,ias,i)=bxcmt(:,:,ias,i)
  enddo
  t1=cb*(bfcmt(3,ia,is)+bfieldc(3))
  bmt(1,:,ias,ndmag)=bmt(1,:,ias,ndmag)+t1/y00
enddo

sv_ubu=0.d0
do ias=1,natmtot
  ic=ias2ic(ias)
  is=ias2is(ias)
! compute radial integrals <u_{l1,io1} | b_{lm3} |  u_{l2,io2}>
  do lm3=0,lmmaxvr
    j1=0
    do l1=0,lmaxvr
      do io1=1,nufr(l1,is) 
        j1=j1+1
        j2=0
        do l2=0,lmaxvr
          do io2=1,nufr(l2,is)
            j2=j2+1
            do i=1,ndmag
              t1=0.d0
              do ir=1,nrmt(is)
                t1=t1+ufr(ir,l1,io1,ic)*ufr(ir,l2,io2,ic)*&
                  bmt(lm3,ir,ias,i)*mt_rw(ir,is)
              enddo
              sv_ubu(lm3,j1,j2,ias,i)=t1
            enddo
          enddo
        enddo
      enddo
    enddo
  enddo
enddo

deallocate(bmt)
return
end subroutine
