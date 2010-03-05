subroutine gen_wann_func(vtrl,ngknr,vgkcnr,wanmt,wanir)
use modmain
implicit none
integer, intent(in) :: vtrl(3)
integer, intent(in) :: ngknr(nkptnrloc)
real(8), intent(in) :: vgkcnr(3,ngkmax,nkptnrloc)
complex(8), intent(out) :: wanmt(lmmaxvr,nrmtmax,natmtot,nspinor,nwann)
complex(8), intent(out) :: wanir(ngrtot,nspinor,nwann)
integer ia,is,ias,ir,ir0,i1,i2,i3,ig,ikloc,ik
integer io,lm,n,ispn,itmp(3)
real(8) v2(3),v3(3),r0,vr0(3),vtrc(3)
complex(8) expikr
logical, external :: vrinmt
wanmt=zzero
wanir=zzero
vtrc(:)=vtrl(1)*avec(:,1)+vtrl(2)*avec(:,2)+vtrl(3)*avec(:,3)
! muffin-tin part
call timer_start(1,reset=.true.)
do ias=1,natmtot
  is=ias2is(ias)
  do ikloc=1,nkptnrloc
    ik=mpi_grid_map(nkptnr,dim_k,loc=ikloc)
    expikr=exp(zi*dot_product(vkcnr(:,ik),vtrc(:)))/nkptnr
    do ir=1,nrmt(is)
      do lm=1,lmmaxvr
        do io=1,nrfl(lm2l(lm),is)
          do ispn=1,nspinor
            do n=1,nwann
              wanmt(lm,ir,ias,ispn,n)=wanmt(lm,ir,ias,ispn,n)+&
                expikr*wann_unkmt(lm,io,ias,ispn,n,ikloc)*&
                urf(ir,lm2l(lm),io,ias)
            enddo !n
          enddo !ispn
        enddo !io
      enddo !lm
    enddo !ir
  enddo !ikloc
enddo !ias 
call timer_stop(1)
! interstitial part
call timer_start(2,reset=.true.)
ir=0
do i3=0,ngrid(3)-1
  v2(3)=dble(i3)/dble(ngrid(3))
  do i2=0,ngrid(2)-1
    v2(2)=dble(i2)/dble(ngrid(2))
    do i1=0,ngrid(1)-1
      v2(1)=dble(i1)/dble(ngrid(1))
      ir=ir+1
      call r3mv(avec,v2,v3)
      if (.not.vrinmt(v3,is,ia,itmp,vr0,ir0,r0)) then
        v3(:)=v3(:)+vtrc(:)
        do ikloc=1,nkptnrloc
          do ig=1,ngknr(ikloc)
            expikr=exp(zi*dot_product(vgkcnr(:,ig,ikloc),v3(:)))
            do ispn=1,nspinor
              do n=1,nwann
                wanir(ir,ispn,n)=wanir(ir,ispn,n)+&
                  expikr*wann_unkit(ig,ispn,n,ikloc)
              enddo !in
            enddo !ispn
          enddo !ig
        enddo !ikloc
      endif
    enddo !i1
  enddo !i2
enddo !i3
wanir(:,:,:)=wanir(:,:,:)/sqrt(omega)/nkptnr
call timer_stop(2)
call mpi_grid_reduce(wanmt(1,1,1,1,1),lmmaxvr*nrmtmax*natmtot*nspinor*nwann,&
  dims=(/dim_k/))
call mpi_grid_reduce(wanir(1,1,1),ngrtot*nspinor*nwann,dims=(/dim_k/))    
return
end

