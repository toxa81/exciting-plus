subroutine seceqnsv_init
use modmain
use mod_sic
implicit none
real(8), allocatable :: bmt(:,:,:,:)
real(8) cb,t1
integer i,ias,ia,is,ic,l1,l2,l3,m1,m2,m3,lm3,lm1,lm2,io1,io2
integer ir,j1,j2
real(8) r2(nrmtmax),fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)
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
  do lm3=1,lmmaxvr
    j1=0
    do l1=0,lmaxapw
      do io1=1,nufr(l1,is) 
        j1=j1+1
        j2=0
        do l2=0,lmaxapw
          do io2=1,nufr(l2,is)
            j2=j2+1
            do i=1,ndmag
              do ir=1,nrmt(is)
                fr(ir)=ufr(ir,l1,io1,ic)*ufr(ir,l2,io2,ic)*&
                  bmt(lm3,ir,ias,i)*(spr(ir,is)**2)
              enddo
              call fderiv(-1,nrmt(is),spr(:,is),fr,gr,cf)
              sv_ubu(lm3,j1,j2,ias,i)=gr(nrmt(is))
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
